import iTrace;
using namespace std;
using namespace numbers;

// Simulate ender eye measurement, with exact error distribution.
// First, strongholds are generated with the preset world seed.
// Next, a random position and measurement error are picked.
// Then, the nearest stronghold is found and error is applied.
// Finally, the command for iTrace is generated to execute.
// I've done massive tests by this, which is impossible by myself.
int main() {
	constexpr long long Base = MC_1_16, Seed = -1236314517;
	constexpr double Emean = 0, Esigma = 0.004, Count = 32;
	Generator Source; StrongholdIter Target;
	setupGenerator(&Source, Base, false);
	applySeed(&Source, DIM_OVERWORLD, Seed);
	initFirstStronghold(&Target, Base, Seed);
	vector<pair<double, double>> data;
	double Offset = Base < MC_1_19 ? Base < MC_1_8 ? 0 : 4 : -4;
	while (nextStronghold(&Target, &Source) > 0)
		data.emplace_back(Target.pos.x + Offset, Target.pos.z + Offset);
	iTrace Instance; string Input;
	ofstream save("data.txt", ios::app);
	Instance(format("VER {0}", mc2str(Base)));
	Instance(format("ERR {0} {1}", Emean, Esigma));
	save << Instance("CHECK") << endl;
	default_random_engine RNG(Seed);
	for (double Index = 0; Index < Count; Index++) {
		double PosX = uniform_int_distribution(-25000, 25000)(RNG);
		double PosZ = uniform_int_distribution(-25000, 25000)(RNG);
		double Error = normal_distribution(Emean, Esigma)(RNG);
		double Yaw = 0, Dmin = +numeric_limits<double>::infinity();
		for (const auto& str : data) {
			double Dist = hypot(str.first - PosX, str.second - PosZ);
			double Angle = 180/pi * atan2(str.second - PosZ, str.first - PosX) - 90;
			if (Dist < Dmin) Dmin = Dist, Yaw = remainder(Angle + Error, 360);
		}
		if (Base < MC_1_13) Input = format("ADD {0:.2f} {1:.2f} {2:.1f}", PosX, PosZ, Yaw);
		else Input = format("/execute in minecraft:overworld run tp @s {0:.2f} 240.00 {1:.2f} {2:.2f} -32.00", PosX, PosZ, Yaw);
		save << Input << endl << Instance(Input) << endl << "CLEAR" << endl << Instance("CLEAR") << endl;
	}
}