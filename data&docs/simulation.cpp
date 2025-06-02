import iTrace;
import "cubiomes/finders.h";
using namespace std;
using namespace numbers;

// Simulate ender eye measurement, with exact error distribution.
// First, strongholds are generated with the world seed.
// Then, a random position and measurement error are picked.
// Finally, the command for iTrace is generated to execute.
// Only MC 1.9+ is supported intentionally, as is Ninjabrain Bot.
// I've taken massive tests by this, which is impossible by myself.
class Simulator {

private:

	struct { double PosX, PosZ; } Map[128];

	double Esum, Count, PosMid;

public:

	Simulator(MCVersion Base, long long Seed) {
		static Generator World;
		static StrongholdIter Target;
		Esum = 0, Count = 0;
		PosMid = Base < MC_1_19 ? 4 : -4;
		setupGenerator(&World, Base, false);
		applySeed(&World, DIM_OVERWORLD, Seed);
		initFirstStronghold(&Target, Base, Seed);
		for (auto& str : Map) {
			nextStronghold(&Target, &World);
			str.PosX = Target.pos.x;
			str.PosZ = Target.pos.z;
		}
	}

	string solve(double Emean, double Evar, bool useDebug) {
		static default_random_engine RNG(time(nullptr));
		double Error = normal_distribution(Emean, Evar)(RNG);
		double Angle = uniform_real_distribution(-pi, pi)(RNG);
		double Radius = uniform_int_distribution(0, 24000)(RNG);
		double PosX = 0.01 * round(100 * Radius * cos(Angle)) - PosMid;
		double PosZ = 0.01 * round(100 * Radius * sin(Angle)) - PosMid;
		double Yaw = 0, Dmin = numeric_limits<double>::max();
		for (const auto& str : Map) {
			double Dist = hypot(PosX - str.PosX, PosZ - str.PosZ);
			double Angle = 180/pi * atan2(str.PosZ - PosZ, str.PosX - PosX) - 90;
			if (Dist < Dmin)
				Dmin = Dist, Yaw = remainder(Angle + Error, 360);
		}
		if (useDebug) return format("ADD {0:.2f} {1:.2f} {2:.1f}", PosX + PosMid, PosZ + PosMid, Yaw);
		else return format("/execute in minecraft:overworld run tp @s {0:.2f} 256.00 {1:.2f} {2:.2f} -32.00", PosX + PosMid, PosZ + PosMid, Yaw);
	}

	string calib(double Emean, double Evar, bool useDebug) {
		static default_random_engine RNG(time(nullptr));
		const auto& Target = Map[0];
		double Error = normal_distribution(Emean, Evar)(RNG);
		double Angle = uniform_real_distribution(-pi, pi)(RNG);
		double PosX = 0.01 * round(100 * (Target.PosX + 64 * cos(Angle))) - PosMid;
		double PosZ = 0.01 * round(100 * (Target.PosZ + 64 * sin(Angle))) - PosMid;
		if (useDebug) {
			double Angle = 180/pi * atan2(Target.PosZ - PosZ, Target.PosX - PosX) - 90;
			double Yaw = 0.1 * round(10 * remainder(Angle + Error, 360));
			double Error = remainder(Yaw - Angle, 360);
			Esum += Error * Error, Count++;
			return format("ADD {0:.2f} {1:.2f} {2:.1f}", PosX + PosMid, PosZ + PosMid, Yaw);
		}
		else {
			double Angle = 180/pi * atan2(Target.PosZ - PosZ, Target.PosX - PosX) - 90;
			double Yaw = 0.01 * round(100 * remainder(Angle + Error, 360));
			double Error = remainder(Yaw - Angle, 360);
			Esum += Error * Error, Count++;
			return format("/execute in minecraft:overworld run tp @s {0:.2f} 256.00 {1:.2f} {2:.2f} -32.00", PosX + PosMid, PosZ + PosMid, Yaw);
		}
	}

	double NBcalib() const { return sqrt(Esum / (Count - 1)); }

};

int main() {
	iTrace Instance; Simulator World(MC_1_16, -1236314517);
	auto execute = [&Instance](const string& Input) { cout << Input << endl << Instance(Input) << endl; };
	execute("VER 1.16"), execute("CAL -1236314517");
	for (int i = 0; i < 32; i++)
		execute(World.calib(0, 0.004, false));
	execute("CAL 0"), execute("CLEAR"), execute("ERR 0 0.004");
	cout << format("NB Standard Deviation: {0:.4f}", World.NBcalib()) << endl;
	for (int i = 0; i < 32; i++)
		execute(World.solve(0, 0.004, false)), execute("CLEAR");
	system("PAUSE");
}