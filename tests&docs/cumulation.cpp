import iTrace;
using namespace std;

// Measure the cumulative distribution of stronghold, with 100M samples.
// Only take radius as the factor, while angle is discarded for simplicity.
// Each ring shares the same distribution, as proven in previous density tests.
// Also versions have been clustered into group 1.0-1.8, 1.9-1.12 and 1.13-since.
// The measurement data are directly used in iTrace in raw. I think it works well.
// Tested version 1.7.10, 1.12.2 and 1.16.5, which represents each group.
int main() {
	constexpr long long Base = MC_1_16, Width = 1, Rmax = 2880, Smax = 1E8;
	Generator Source; StrongholdIter Target;
	long long Progress = 0, Count[Rmax / Width]{};
	for (long long Seed = 1; Seed <= Smax; Seed++) {
		if (Seed == Progress * Smax / 100) cout << '|', Progress++;
		setupGenerator(&Source, Base, false);
		applySeed(&Source, DIM_OVERWORLD, Seed);
		initFirstStronghold(&Target, Base, Seed);
		nextStronghold(&Target, &Source);
		long long Radius = hypot(Target.pos.x + 4, Target.pos.z + 4);
		for (size_t Index = Radius / Width; Index < Rmax / Width; Count[Index++]++);
	}
	ofstream save("data.csv", ios::app);
	save << mc2str(Base) << ',' << Smax << endl;
	for (size_t Index = 0; Index < Rmax / Width; Index++)
		save << Index * Width << ',' << Count[Index] << endl;
}