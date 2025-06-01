import std;
import "cubiomes/finders.h";
using namespace std;

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