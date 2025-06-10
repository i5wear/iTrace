import std;
import "cubiomes/finders.h";
using namespace std;

// Measure the cumulative distribution of stronghold, with 100M samples.
// Only take radius as the factor, while angle is discarded for simplicity.
// Each ring shares the same distribution, as proven in previous density tests.
// Also versions have been clustered into group 1.0-1.8, 1.9-1.12 and 1.13-since.
// The measurement data are directly used in iTrace in raw. I think it works well.
// Tested version 1.7.10, 1.12.2 and 1.16.5, which represents each group.
int main() {
	constexpr long long Base = MC_1_16, Rmin = 1200, Rmax = 2880, Width = 1, Smin = 1, Smax = 1E8, Trial = 1;
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
			for (long long N = (Radius - Rmin) / Width; N < (Rmax - Rmin) / Width; N++)
				Function[N]++;
		}
	}
	string Output = "Radius,Samples\n";
	for (long long N = 0; N < (Rmax - Rmin) / Width; N++)
		Output += format("{0},{1}\n", Rmin + N * Width, Function[N]);
	ofstream("data.csv") << Output;
}