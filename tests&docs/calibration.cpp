import iTrace;
using namespace std;
using namespace numbers;

// Simulate error calibration, with exact error distribution.
// First, strongholds are generated with the preset world seed.
// Next, a random error is picked and position is given by iTrace.
// Then, the nearest stronghold is found and error is applied.
// Finally, the command for iTrace is generated to execute.
// Also simulate Ninjabrain Bot behavior to share the same dataset.
int main() {
	constexpr long long Base = MC_1_16, Seed = -1236314517;
	constexpr double Emean = 0, Esigma = 0.004, Count = 32;
	Generator Source; StrongholdIter Target;
	setupGenerator(&Source, Base, false);
	applySeed(&Source, DIM_OVERWORLD, Seed);
	initFirstStronghold(&Target, Base, Seed);
	vector<pair<double, double>> Data;
	double Offset = Base < MC_1_19 ? Base < MC_1_8 ? 0 : 4 : -4;
	while (nextStronghold(&Target, &Source) > 0)
		Data.emplace_back(Target.pos.x + Offset, Target.pos.z + Offset);
	iTrace Instance; string Input;
	ofstream save("data.txt", ios::app);
	Instance(format("VER {0}", mc2str(Base)));
	Instance(format("ERR {0} {1}", Emean, Esigma));
	save << Instance("CHECK") << endl;
	const regex Pattern("#[0-9]+: (\\S+) → ERR: (\\S+) ± (\\S+)\n/tp (\\S+) (\\S+) (\\S+)\n", regex::icase);
	string Output; smatch Value; double Esum2 = 0;
	Input = format("CAL {0}", Seed), Output = Instance(Input);
	regex_match(Output, Value, Pattern);
	save << Input << endl << Output << endl;
	default_random_engine RNG(Seed);
	for (double Index = 0; Index < Count; Index++) {
		double PosX = stod(Value[4]), PosZ = stod(Value[6]);
		double Error = normal_distribution(Emean, Esigma)(RNG);
		double Yaw = 0, Dmin = +numeric_limits<double>::infinity();
		for (const auto& Pair : Data) {
			double Distance = hypot(Pair.first - PosX, Pair.second - PosZ);
			double Angle = 180/pi * atan2(Pair.second - PosZ, Pair.first - PosX) - 90;
			if (Distance < Dmin) Dmin = Distance, Yaw = remainder(Angle + Error, 360);
		}
		if (Base < MC_1_13) Input = format("ADD {0:.2f} {1:.2f} {2:.1f}", PosX, PosZ, Yaw);
		else Input = format("/execute in minecraft:overworld run tp @s {0:.2f} 240.00 {1:.2f} {2:.2f} -32.00", PosX, PosZ, Yaw);
		Output = Instance(Input), regex_match(Output, Value, Pattern), Esum2 += stod(Value[1]) * stod(Value[1]);
		save << Input << endl << Output << endl << format("Ninjabrain Bot: 0.0000 ± {0:.4f}\n", sqrt(Esum2 / Index)) << endl;
	}
}