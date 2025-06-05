import iTrace;
import "cubiomes/finders.h";
import "cubiomes/util.h";
using namespace std;
using namespace numbers;

// Simulate ender eye measurement, with exact error distribution.
// First, strongholds are generated with the preset world seed.
// Then, a random position and measurement error are picked.
// Finally, the command for iTrace is generated to execute.
// Only version 1.9+ is supported intentionally, as Ninjabrain Bot.
// I've taken massive tests by this, which is impossible by myself.
class Simulator {

private:

	struct { double PosX, PosZ; } data[128];

	default_random_engine RNG{ random_device()() };

	double Esum2{ 0 }, Count{ 0 };

public:

	Simulator(int Base, long long Seed) {
		thread_local Generator World;
		thread_local StrongholdIter Target;
		setupGenerator(&World, Base, false);
		applySeed(&World, DIM_OVERWORLD, Seed);
		initFirstStronghold(&Target, Base, Seed);
		double Offset = Base < MC_1_19 ? 4 : -4;
		for (auto& str : data) {
			nextStronghold(&Target, &World);
			str.PosX = Target.pos.x + Offset;
			str.PosZ = Target.pos.z + Offset;
		}
	}

	string solve(double Emean, double Evar, bool Legacy) {
		double Error = normal_distribution(Emean, Evar)(RNG);
		double Angle = uniform_real_distribution(-pi, pi)(RNG);
		double Radius = uniform_int_distribution(0, 24000)(RNG);
		double PosX = 0.01 * round(100 * Radius * cos(Angle));
		double PosZ = 0.01 * round(100 * Radius * sin(Angle));
		double Yaw = 0, Dmin = +numeric_limits<double>::infinity();
		for (const auto& str : data) {
			double Dist = hypot(PosX - str.PosX, PosZ - str.PosZ);
			double Angle = 180/pi * atan2(str.PosZ - PosZ, str.PosX - PosX) - 90;
			if (Dist < Dmin) Dmin = Dist, Yaw = remainder(Angle + Error, 360);
		}
		if (Legacy) return format("ADD {0:.2f} {1:.2f} {2:.1f}", PosX, PosZ, Yaw);
		else return format("/execute in minecraft:overworld run tp @s {0:.2f} 250.00 {1:.2f} {2:.2f} -32.00", PosX, PosZ, Yaw);
	}

	string calib(double Emean, double Evar, bool Legacy) {
		const auto& Target{ data[0] };
		double Error = normal_distribution(Emean, Evar)(RNG);
		double Angle = uniform_real_distribution(-pi, pi)(RNG);
		double PosX = 0.01 * round(100 * (Target.PosX + 60 * cos(Angle)));
		double PosZ = 0.01 * round(100 * (Target.PosZ + 60 * sin(Angle)));
		if (Legacy) {
			double Angle = 180/pi * atan2(Target.PosZ - PosZ, Target.PosX - PosX) - 90;
			double Yaw = 0.1 * round(10 * remainder(Angle + Error, 360));
			double Error = remainder(Yaw - Angle, 360);
			Esum2 += Error * Error, Count++;
			return format("ADD {0:.2f} {1:.2f} {2:.1f}", PosX, PosZ, Yaw);
		}
		else {
			double Angle = 180/pi * atan2(Target.PosZ - PosZ, Target.PosX - PosX) - 90;
			double Yaw = 0.01 * round(100 * remainder(Angle + Error, 360));
			double Error = remainder(Yaw - Angle, 360);
			Esum2 += Error * Error, Count++;
			return format("/execute in minecraft:overworld run tp @s {0:.2f} 250.00 {1:.2f} {2:.2f} -32.00", PosX, PosZ, Yaw);
		}
	}

	double NBcalib() const { return sqrt(Esum2 / (Count - 1)); }

};

int main() {
	constexpr int Base = MC_1_16, Seed = -1236314517, Trial = 30;
	constexpr double Emean = 0, Evar = 0.004;
	iTrace Instance; Simulator World{ Base, Seed };
	auto execute = [&Instance](const string& Input)
		{ cout << Input << endl << Instance(Input) << endl; };
	execute(format("VER {0}", mc2str(Base))), execute(format("CAL {0}", Seed));
	for (int i = 0; i < Trial; i++)
		execute(World.calib(Emean, Evar, false));
	cout << format("NB Standard Deviation: {0:.4f}\n", World.NBcalib()) << endl;
	execute("CAL 0"), execute(format("ERR {0} {1}", Emean, Evar));
	for (int i = 0; i < Trial; i++)
		execute("CLEAR"), execute(World.solve(Emean, Evar, false));
	system("PAUSE");
}