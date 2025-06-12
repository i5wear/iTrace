import iTrace;
import "cubiomes/finders.h";
import "cubiomes/util.h";
using namespace std;
using namespace numbers;
#define VALUE " ([+-]?[0-9]+(?:[.][0-9]+)?)"

// Simulate error calibration, with exact error distribution.
// First, strongholds are generated with the preset world seed.
// Next, a random error is picked and position is given by iTrace.
// Then, the nearest stronghold is found and error is applied.
// Finally, the command for iTrace is generated to execute.
// Also simulate ninjabrain bot behavior to share the same dataset.
int main() {
	constexpr long long Base = MC_1_16, Seed = -1236314517;
	constexpr double Emean = 0, Evar = 0.004, Count = 32;
	default_random_engine RNG(Seed);
	Generator Source; StrongholdIter Target;
	setupGenerator(&Source, Base, false);
	applySeed(&Source, DIM_OVERWORLD, Seed);
	initFirstStronghold(&Target, Base, Seed);
	vector<pair<double, double>> data;
	double Esum2 = 0, Offset = Base < MC_1_19 ? Base < MC_1_8 ? 0 : 4 : -4;
	while (nextStronghold(&Target, &Source) > 0)
		data.emplace_back(Target.pos.x + Offset, Target.pos.z + Offset);
	regex Feedback("#[0-9]+:" VALUE ", MEAN:" VALUE ", SD:" VALUE, regex::icase);
	regex Command("/tp" VALUE VALUE VALUE, regex::icase);
	iTrace Instance; smatch Value;
	ofstream save("data.txt", ios::noreplace);
	Instance(format("VER {0}", mc2str(Base)));
	Instance(format("ERR {0} {1}", Emean, Evar));
	save << Instance("CHECK") << endl;
	string Input = format("CAL {0}", Seed);
	string Output = Instance(Input);
	save << Input << endl << Output << endl;
	for (double Index = 0; Index < Count; Index++) {
		regex_search(Output, Value, Command);
		double PosX = stod(Value[1]), PosZ = stod(Value[3]);
		double Error = normal_distribution(Emean, Evar)(RNG);
		double Yaw = 0, Dmin = +numeric_limits<double>::infinity();
		for (const auto& str : data) {
			double Dist = hypot(str.first - PosX, str.second - PosZ);
			double Angle = 180/pi * atan2(str.second - PosZ, str.first - PosX) - 90;
			if (Dist < Dmin) Dmin = Dist, Yaw = remainder(Angle + Error, 360);
		}
		if (Base < MC_1_13) Input = format("ADD {0:.2f} {1:.2f} {2:.1f}", PosX, PosZ, Yaw);
		else Input = format("/execute in minecraft:overworld run tp @s {0:.2f} 240.00 {1:.2f} {2:.2f} -32.00", PosX, PosZ, Yaw);
		Output = Instance(Input);
		regex_search(Output, Value, Feedback);
		Esum2 += stod(Value[1]) * stod(Value[1]);
		save << Input << endl << Output << endl << format("SD(Ninjabrain Bot): {0:.4f}\n", sqrt(Esum2 / Index)) << endl;
	}
}