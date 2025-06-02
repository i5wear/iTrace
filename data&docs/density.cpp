import std;
import "cubiomes/finders.h";
using namespace std;

// Measure the probability density of stronghold, with 30M samples.
// Since calculation for MC 1.18+ is very slow, it's reduced to 3M samples.
// Only take radius as the factor, while angle is discarded for simplicity.
// Cluster MC versions into different groups by the measurement data.
// Also take 'global' tests, to prove each ring shares the same distribution.
// Tested version 1.0, 1.9, 1.11, 1.12, 1.13, 1.15, 1.17 and 1.20.
int main() {
	constexpr int Base = MC_1_20, Rmin = 1200, Rmax = 2880, Width = 16, Smin = 1, Smax = 1E6, Trial = 3;
	int Function[(Rmax - Rmin) / Width] = { 0 };
	Generator World;
	StrongholdIter Target;
	setupGenerator(&World, Base, false);
	for (int Seed = Smin, Progress = 0; Seed <= Smax; Seed++) {
		if (Seed - Smin == Progress * (Smax - Smin) / 100) {
			system("CLS");
			cout << Progress << "%\n";
			Progress++;
		}
		initFirstStronghold(&Target, Base, Seed);
		applySeed(&World, DIM_OVERWORLD, Seed);
		for (int i = 0; i < Trial; i++) {
			nextStronghold(&Target, &World);
			Function[int(hypot(Target.pos.x, Target.pos.z) - Rmin) / Width]++;
		}
	}
	string Output("Range,Samples\n");
	for (int N = 0; N < (Rmax - Rmin) / Width; N++)
		Output.append(format("{}, {}\n", Rmin + N * Width, Function[N]));
	ofstream("data.csv") << Output;
	system("PAUSE");
}