import iTrace;
using namespace std;
using namespace numbers;

/* Release 1 - iTrace Console User Interface */
int main() {
    iTrace Instance; string Input;
    ofstream save("history.log", ios::app);
    for (ifstream load("history.log"); getline(load, Input).good(); Instance(Input));
    cout << Instance("CHECK") << endl;
    while (getline(cin, Input).good() and not Input.empty()) {
        try { cout << Instance(Input) << endl; }
        catch (...) { continue; }
        save << Input << endl;
    }
} /**/

/* Test 1 - Density of Stronghold Distribution /
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
        save << (Index + 1) * Width << ',' << Count[Index] << endl;
} /**/

/* Test 2 - Cumulative Stronghold Distribution /
int main() {
    constexpr long long Base = MC_1_16, Width = 1, Rmax = 2880, Smax = 1E8;
    Generator Source; StrongholdIter Target;
    long long Progress = 1, Count[Rmax / Width]{};
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
        save << (Index + 1) * Width << ',' << Count[Index] << endl;
} /**/

/* Test 3 - Simulation of Endereye Measurement /
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
    default_random_engine RNG(Seed);
    for (double Index = 0; Index < Count; Index++) {
        double PosX = uniform_int_distribution(-25000, 25000)(RNG);
        double PosZ = uniform_int_distribution(-25000, 25000)(RNG);
        double Error = normal_distribution(Emean, Esigma)(RNG);
        double Yaw = 0, Dmin = +numeric_limits<double>::infinity();
        for (const auto& Pair : Data) {
            double Distance = hypot(Pair.first - PosX, Pair.second - PosZ);
            double Angle = 180 / pi * atan2(Pair.second - PosZ, Pair.first - PosX) - 90;
            if (Distance < Dmin) Dmin = Distance, Yaw = remainder(Angle + Error, 360);
        }
        if (Base < MC_1_13) Input = format("ADD {0:.2f} {1:.2f} {2:.1f}", PosX, PosZ, Yaw);
        else Input = format("/execute in minecraft:overworld run tp @s {0:.2f} 240.00 {1:.2f} {2:.2f} -32.00", PosX, PosZ, Yaw);
        save << Input << endl << Instance(Input) << endl << "CLEAR" << endl << Instance("CLEAR") << endl;
    }
} /**/

/* Test 4 - Simulation of Error Calibration /
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
            double Angle = 180 / pi * atan2(Pair.second - PosZ, Pair.first - PosX) - 90;
            if (Distance < Dmin) Dmin = Distance, Yaw = remainder(Angle + Error, 360);
        }
        if (Base < MC_1_13) Input = format("ADD {0:.2f} {1:.2f} {2:.1f}", PosX, PosZ, Yaw);
        else Input = format("/execute in minecraft:overworld run tp @s {0:.2f} 240.00 {1:.2f} {2:.2f} -32.00", PosX, PosZ, Yaw);
        Output = Instance(Input), regex_match(Output, Value, Pattern), Esum2 += stod(Value[1]) * stod(Value[1]);
        save << Input << endl << Output << endl << format("Ninjabrain Bot: 0.0000 ± {0:.4f}\n", sqrt(Esum2 / Index)) << endl;
    }
} /**/