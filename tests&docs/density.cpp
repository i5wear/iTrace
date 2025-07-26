import iTrace;
using namespace std;

// Measure the probability density of stronghold, with 12.8M samples.
// Only take radius as the factor, while angle is discarded for simplicity.
// Cluster game versions into different groups by the measurement data.
// Also aim to prove that each ring shares the same distribution.
// Tested final versions of 1.0, 1.8, 1.9, 1.12, 1.13, 1.17 and 1.20.
int main() {
	constexpr long long Base = MC_1_20, Width = 128, Rmax = 25600, Smax = 1E5;
	Generator Source; StrongholdIter Target;
	long long Progress = 1, Count[Rmax / Width]{};
	for (long long Seed = 1; Seed <= Smax; Seed++) {
		if (Seed == Progress * Smax / 100) cout << '|', Progress++;
		setupGenerator(&Source, Base, false);
		applySeed(&Source, DIM_OVERWORLD, Seed);
		initFirstStronghold(&Target, Base, Seed);
		while (nextStronghold(&Target, &Source) > 0)
			Count[size_t(hypot(Target.pos.x + 4, Target.pos.z + 4) / Width)]++;
	}
	ofstream save("data.csv", ios::app);
	save << mc2str(Base) << ',' << Smax << endl;
	for (size_t Index = 0; Index < Rmax / Width; Index++)
		save << Index * Width << ',' << Count[Index] << endl;
}