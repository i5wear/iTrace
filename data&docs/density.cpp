import std;
import "cubiomes/finders.h";
using namespace std;

// Measure the probability density of stronghold, with 12.8M samples.
// Since calculation on MC 1.18+ is very slow, it's reduced to 128K.
// Only take radius as the factor, while angle is discarded for simplicity.
// Cluster game versions into different groups by the measurement data.
// Also aim to prove that each ring shares the same distribution.
// Tested final versions of 1.0, 1.8, 1.9, 1.12, 1.13, 1.17, 1.18 and 1.20.
int main() {
	constexpr long long Base = MC_1_20, Rmin = 0, Rmax = 25600, Width = 128, Smin = 1, Smax = 1E3, Trial = 128;
	Generator World; StrongholdIter Target;
	long long Progress = 0, Function[(Rmax - Rmin) / Width] = {0};
	for (long long Seed = Smin; Seed <= Smax; Seed++) {
		if (Seed - Smin == Progress * (Smax - Smin) / 100)
			system("CLS"), cout << format("{0}%\n", ++Progress);
		setupGenerator(&World, Base, false);
		applySeed(&World, DIM_OVERWORLD, Seed);
		initFirstStronghold(&Target, Base, Seed);
		for (long long i = 0; i < Trial; i++) {
			nextStronghold(&Target, &World);
			long long Radius = hypot(Target.pos.x + 4, Target.pos.z + 4);
			Function[(Radius - Rmin) / Width]++;
		}
	}
	string Output = "Radius,Samples\n";
	for (long long N = 0; N < (Rmax - Rmin) / Width; N++)
		Output += format("{0},{1}\n", Rmin + N * Width, Function[N]);
	ofstream("data.csv") << Output;
}