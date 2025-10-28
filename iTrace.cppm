export module iTrace;
export import "cubiomes/finders.h";
export import "cubiomes/util.h";
export import std;
using namespace std;
using namespace numbers;

export class iTrace {

protected:

    struct Constants {
        struct Ring { double Rmin, Rmax, Count; const double* Distr; };
        double Chunk, PosLoc, PosMid; vector<Ring> Data;
        static const double TABLE7[], TABLE12[], TABLE16[];
        Constants(long long Base) {
            if (Base < MC_1_9) {
                Chunk = 16, PosLoc = 4, PosMid = 4;
                Data.emplace_back(460, 1332, 3, TABLE7);
            }
            else if (Base < MC_1_13) {
                Chunk = 16, PosLoc = 4, PosMid = 0;
                Data = {
                    {  1228,  2868,  3, TABLE12 },
                    {  4300,  5940,  6, TABLE12 },
                    {  7372,  9012, 10, TABLE12 },
                    { 10444, 12084, 15, TABLE12 },
                    { 13516, 15156, 21, TABLE12 },
                    { 16588, 18228, 28, TABLE12 },
                    { 19660, 21300, 36, TABLE12 },
                    { 22732, 24372,  9, TABLE12 }
                };
            }
            else if (Base < MC_1_19) {
                Chunk = 16, PosLoc = 4, PosMid = 0;
                Data = {
                    {  1228,  2868,  3, TABLE16 },
                    {  4300,  5940,  6, TABLE16 },
                    {  7372,  9012, 10, TABLE16 },
                    { 10444, 12084, 15, TABLE16 },
                    { 13516, 15156, 21, TABLE16 },
                    { 16588, 18228, 28, TABLE16 },
                    { 19660, 21300, 36, TABLE16 },
                    { 22732, 24372,  9, TABLE16 }
                };
            }
            else {
                Chunk = 16, PosLoc = 4, PosMid = 8;
                Data = {
                    {  1228,  2868,  3, TABLE16 },
                    {  4300,  5940,  6, TABLE16 },
                    {  7372,  9012, 10, TABLE16 },
                    { 10444, 12084, 15, TABLE16 },
                    { 13516, 15156, 21, TABLE16 },
                    { 16588, 18228, 28, TABLE16 },
                    { 19660, 21300, 36, TABLE16 },
                    { 22732, 24372,  9, TABLE16 }
                };
            }
        }
    };

    struct Endereyes {
        struct Line { double PosX, PosZ, Yaw, Range; };
        mutable double Emean, Esigma; vector<Line> Data;
        double calib(const Constants& Base, double PosX, double PosZ) const {
            PosX = Base.Chunk * (floor(PosX / Base.Chunk) + 0.5) - Base.PosMid;
            PosZ = Base.Chunk * (floor(PosZ / Base.Chunk) + 0.5) - Base.PosMid;
            double Error = 0, Esum1 = 0, Esum2 = 0, Rsum2 = 0;
            for (const auto& Line : Data) {
                Error = remainder(Line.Yaw - atan2(PosZ - Line.PosZ, PosX - Line.PosX), 2 * pi);
                Esum1 += Error, Esum2 += Error * Error, Rsum2 += Line.Range * Line.Range;
            }
            if (not Data.empty()) {
                Esum1 /= Data.size(), Esum2 /= Data.size() - 1, Rsum2 /= Data.size();
                Emean = Esum1, Esigma = sqrt(fdim(Esum2, Esum1 * Esum1 + Rsum2 / 3));
            }
            return Error;
        }
        double solve(const Constants& Base, double PosX, double PosZ) const {
            PosX = Base.Chunk * (floor(PosX / Base.Chunk) + 0.5) - Base.PosMid;
            PosZ = Base.Chunk * (floor(PosZ / Base.Chunk) + 0.5) - Base.PosMid;
            double Prob = 0, Radius = hypot(PosX + Base.PosMid, PosZ + Base.PosMid);
            for (const auto& Ring : Base.Data)
                if (Radius > Ring.Rmin and Radius < Ring.Rmax)
                    Prob = Ring.Count * (Ring.Rmax - Ring.Rmin) / (Base.Chunk * Radius);
            double Dmin = +numeric_limits<double>::infinity();
            double Dmax = -numeric_limits<double>::infinity();
            for (const auto& Line : Data) {
                double Error = remainder(Line.Yaw - atan2(PosZ - Line.PosZ, PosX - Line.PosX), 2 * pi);
                double Emin = (Error - Emean - Line.Range) / (sqrt2 * Esigma);
                double Emax = (Error - Emean + Line.Range) / (sqrt2 * Esigma);
                Prob *= Esigma ? Esigma * (erf(Emax) - erf(Emin)) / Line.Range : Emin < 0 and Emax > 0;
                Dmin = fmin(Dmin, hypot(Line.PosX + Base.PosMid, Line.PosZ + Base.PosMid) - hypot(PosX - Line.PosX, PosZ - Line.PosZ));
                Dmax = fmax(Dmax, hypot(Line.PosX + Base.PosMid, Line.PosZ + Base.PosMid) + hypot(PosX - Line.PosX, PosZ - Line.PosZ));
            }
            vector<pair<double, double>> Interval;
            for (const auto& Ring : Base.Data) {
                auto Distr = [&Ring](double Radius) { return Radius < Ring.Rmax ? Radius < Ring.Rmin ? 0 : Ring.Distr[size_t(Radius - Ring.Rmin)] : 1; };
                if (Prob <= 0) break;
                else if (Radius > Ring.Rmin and Radius < Ring.Rmax) {
                    Prob *= Distr(Radius + 0.5 * Base.Chunk) - Distr(Radius - 0.5 * Base.Chunk);
                    for (double Index = 1; Index < Ring.Count; Index++) {
                        double Angle = 2 * pi * Index / Ring.Count + atan2(PosZ, PosX);
                        for (const auto& Line : Data) {
                            double Scale = (PosX + Base.PosMid) * (Line.PosX + Base.PosMid) + (PosZ + Base.PosMid) * (Line.PosZ + Base.PosMid);
                            double Coeff = (Line.PosX + Base.PosMid) * cos(Angle) + (Line.PosZ + Base.PosMid) * sin(Angle);
                            double Delta = Radius * Radius + Coeff * Coeff - 2 * Scale;
                            if (Delta > 0) Interval.emplace_back(Coeff - sqrt(Delta), Coeff + sqrt(Delta));
                        }
                        double Pmult = 1, Rmax = -numeric_limits<double>::infinity();
                        ranges::sort(Interval, ranges::less());
                        for (const auto& Pair : Interval) {
                            Pmult -= Distr(fmax(Rmax, Pair.second)) - Distr(fmax(Rmax, Pair.first));
                            Rmax = fmax(Rmax, Pair.second);
                        }
                        Prob *= Pmult, Interval.clear();
                    }
                }
                else if (Dmax > Ring.Rmin and Dmin < Ring.Rmax) {
                    double Psum = 0, Count = round(pi * (Ring.Rmin + Ring.Rmax) / (Base.Chunk * Ring.Count));
                    for (double Step = 0.5; Step < Count; Step++) {
                        double Prob = 1 / Count;
                        for (double Index = 0; Index < Ring.Count; Index++) {
                            double Angle = 2 * pi * (Index + Step / Count) / Ring.Count;
                            for (const auto& Line : Data) {
                                double Scale = (PosX + Base.PosMid) * (Line.PosX + Base.PosMid) + (PosZ + Base.PosMid) * (Line.PosZ + Base.PosMid);
                                double Coeff = (Line.PosX + Base.PosMid) * cos(Angle) + (Line.PosZ + Base.PosMid) * sin(Angle);
                                double Delta = Radius * Radius + Coeff * Coeff - 2 * Scale;
                                if (Delta > 0) Interval.emplace_back(Coeff - sqrt(Delta), Coeff + sqrt(Delta));
                            }
                            double Pmult = 1, Rmax = -numeric_limits<double>::infinity();
                            ranges::sort(Interval, ranges::less());
                            for (const auto& Pair : Interval) {
                                Pmult -= Distr(fmax(Rmax, Pair.second)) - Distr(fmax(Rmax, Pair.first));
                                Rmax = fmax(Rmax, Pair.second);
                            }
                            Prob *= Pmult, Interval.clear();
                        }
                        Psum += Prob;
                    }
                    Prob *= Psum;
                }
            }
            return Prob;
        }
    };

    struct Stronghold {
        struct Point { double PosX, PosZ, Prob; };
        double Xmean, Xsigma, Zmean, Zsigma; vector<Point> Data;
        Stronghold(long long Base, long long Seed) {
            thread_local Generator Source;
            thread_local StrongholdIter Target;
            setupGenerator(&Source, Base, false);
            applySeed(&Source, DIM_OVERWORLD, Seed);
            initFirstStronghold(&Target, Base, Seed);
            nextStronghold(&Target, &Source);
            Xmean = Target.pos.x, Xsigma = 0;
            Zmean = Target.pos.z, Zsigma = 0;
            Data.emplace_back(Target.pos.x, Target.pos.z, 1);
        }
        Stronghold(const Constants& Base, const Endereyes& Source) {
            vector<size_t> Step; vector<vector<pair<double, double>>> Dataset;
            for (const auto& Line : Source.Data) {
                Step.emplace_back(), Dataset.emplace_back();
                double Radius = hypot(Line.PosX + Base.PosMid, Line.PosZ + Base.PosMid);
                double Dmin = +numeric_limits<double>::infinity();
                double Dmax = +numeric_limits<double>::infinity();
                for (const auto& Ring : Base.Data) {
                    double Angle = pi / Ring.Count;
                    Dmin = fmin(Dmin, Radius < Ring.Rmax ? Radius < Ring.Rmin ? Ring.Rmin - Radius : 0 : Radius - Ring.Rmax);
                    Dmax = fmin(Dmax, hypot(Radius * sin(Angle), fmax(Radius * cos(Angle) - Ring.Rmin, Ring.Rmax - Radius * cos(Angle))));
                }
                double Angle = Line.Yaw - Source.Emean;
                double Amin = Angle - Line.Range - 4 * Source.Esigma;
                double Amax = Angle + Line.Range + 4 * Source.Esigma;
                double Pos1 = Line.PosX + Dmin * cos(Amin), Pos2 = Line.PosX + Dmax * cos(Amin);
                double Pos3 = Line.PosX + Dmin * cos(Amax), Pos4 = Line.PosX + Dmax * cos(Amax);
                double Xmin = fmin(fmin(Pos1, Pos2), fmin(Pos3, Pos4));
                double Xmax = fmax(fmin(Pos1, Pos2), fmax(Pos3, Pos4));
                Xmin = Base.Chunk * (round((Xmin + Base.PosMid) / Base.Chunk) + 0.5) - Base.PosMid;
                Xmax = Base.Chunk * (round((Xmax + Base.PosMid) / Base.Chunk) + 0.5) - Base.PosMid;
                for (double PosX = Xmin; PosX < Xmax; PosX += Base.Chunk) {
                    double Zmin = +numeric_limits<double>::infinity();
                    double Zmax = -numeric_limits<double>::infinity();
                    if (PosX > fmin(Pos1, Pos2) and PosX < fmax(Pos1, Pos2)) {
                        Zmin = fmin(Zmin, Line.PosZ + (PosX - Line.PosX) * tan(Amin));
                        Zmax = fmax(Zmax, Line.PosZ + (PosX - Line.PosX) * tan(Amin));
                    }
                    if (PosX > fmin(Pos3, Pos4) and PosX < fmax(Pos3, Pos4)) {
                        Zmin = fmin(Zmin, Line.PosZ + (PosX - Line.PosX) * tan(Amax));
                        Zmax = fmax(Zmax, Line.PosZ + (PosX - Line.PosX) * tan(Amax));
                    }
                    if (PosX > fmin(Pos1, Pos3) and PosX < fmax(Pos1, Pos3)) {
                        Zmin = fmin(Zmin, Line.PosZ + Dmin / sin(Angle) - (PosX - Line.PosX) / tan(Angle));
                        Zmax = fmax(Zmax, Line.PosZ + Dmin / sin(Angle) - (PosX - Line.PosX) / tan(Angle));
                    }
                    if (PosX > fmin(Pos2, Pos4) and PosX < fmax(Pos2, Pos4)) {
                        Zmin = fmin(Zmin, Line.PosZ + Dmax / sin(Angle) - (PosX - Line.PosX) / tan(Angle));
                        Zmax = fmax(Zmax, Line.PosZ + Dmax / sin(Angle) - (PosX - Line.PosX) / tan(Angle));
                    }
                    Zmin = Base.Chunk * (round((Zmin + Base.PosMid) / Base.Chunk) + 0.5) - Base.PosMid;
                    Zmax = Base.Chunk * (round((Zmax + Base.PosMid) / Base.Chunk) + 0.5) - Base.PosMid;
                    for (double PosZ = Zmin; PosZ < Zmax; PosZ += Base.Chunk)
                        Dataset.back().emplace_back(PosX, PosZ);
                }
            }
            double Xsum1 = 0, Xsum2 = 0, Zsum1 = 0, Zsum2 = 0, Psum = 0;
            if (not Dataset.empty()) LOOP: {
                size_t IDmin = 0, IDmax = 0;
                for (size_t Index = 0; Index < Dataset.size(); Index++) {
                    if (Step[Index] == Dataset[Index].size()) goto EXIT;
                    else if (Dataset[Index][Step[Index]] < Dataset[IDmin][Step[IDmin]]) Step[Index]++, IDmin = Index;
                    else if (Dataset[Index][Step[Index]] < Dataset[IDmax][Step[IDmax]]) Step[Index]++;
                    else if (Dataset[Index][Step[Index]] > Dataset[IDmax][Step[IDmax]]) Step[IDmax]++, IDmax = Index;
                }
                if (not IDmin and not IDmax) {
                    double PosX = Base.Chunk * floor(Dataset[0][Step[0]]. first / Base.Chunk) + Base.PosLoc;
                    double PosZ = Base.Chunk * floor(Dataset[0][Step[0]].second / Base.Chunk) + Base.PosLoc;
                    double Prob = Source.solve(Base, PosX, PosZ);
                    for (size_t Index = 0; Index < Dataset.size(); Step[Index++]++);
                    if (Prob > 0) {
                        Xsum1 += PosX * Prob, Xsum2 += PosX * PosX * Prob;
                        Zsum1 += PosZ * Prob, Zsum2 += PosZ * PosZ * Prob;
                        Psum += Prob, Data.emplace_back(PosX, PosZ, Prob);
                    }
                }
                goto LOOP;
            } EXIT:
            if (Psum > 0) {
                Xsum1 /= Psum, Xsum2 /= Psum, Zsum1 /= Psum, Zsum2 /= Psum;
                for (volatile auto& Point : Data) Point.Prob /= Psum;
                ranges::sort(Data, ranges::greater(), &Point::Prob);
            }
            Xmean = Xsum1, Xsigma = sqrt(fdim(Xsum2, Xsum1 * Xsum1));
            Zmean = Zsum1, Zsigma = sqrt(fdim(Zsum2, Zsum1 * Zsum1));
        }
    };

private:

    long long Base = MC_NEWEST, Seed = 0;

    Endereyes Source = { pi/720, pi/3600 };

public:

    string operator()(const string& Input) {
        static const regex Pattern[] = {
            regex("CHECK", regex::icase),
            regex("CLEAR", regex::icase),
            regex("VER (\\S+)", regex::icase),
            regex("CAL (\\S+)", regex::icase),
            regex("ERR (\\S+) (\\S+)", regex::icase),
            regex("ADD (\\S+) (\\S+) (\\S+)", regex::icase),
            regex("/execute in minecraft:overworld run tp @s (\\S+) (\\S+) (\\S+) (\\S+) (\\S+)", regex::icase)
        };
        thread_local default_random_engine RNG;
        size_t Index; smatch Value; string Output;
        for (Index = 0; Index < 7; Index++)
            if (regex_match(Input, Value, Pattern[Index])) break;
        switch (Index) {
        case 0: Index = 0; break;
        case 1: Seed = 0, Source.Data.clear(); break;
        case 2: Base = str2mc(Value[1].str().data()); break;
        case 3: Seed = stoll(Value[1]), RNG.seed(stoll(Value[1])); break;
        case 4: Source.Emean = pi/180 * stod(Value[1]), Source.Esigma = pi/180 * stod(Value[2]); break;
        case 5: Source.Data.emplace_back(stod(Value[1]), stod(Value[2]), pi/180 * (stod(Value[3]) + 90), pi/3600); break;
        case 6: Source.Data.emplace_back(stod(Value[1]), stod(Value[3]), pi/180 * (stod(Value[4]) + 90), pi/36000); break;
        default: throw exception("invalid iTrace input");
        }
        if (not Index) {
            Output += format("VER: {0}  SEED: {1}  ERR: {2:.4f} ± {3:.4f}\n", mc2str(Base), Seed, 180/pi * Source.Emean, 180/pi * Source.Esigma);
            for (const auto& Line : Source.Data)
                Output += format("#{0}: ({1:.1f}, {2:.1f}, {3:.2f} ± {4:.1g})\n", ++Index, Line.PosX, Line.PosZ, 180/pi * Line.Yaw - 90, 180/pi * Line.Range);
        }
        else if (Seed) {
            const auto& Target = Stronghold(Base, Seed).Data.back();
            double Angle = uniform_real_distribution(-pi, pi)(RNG);
            double Error = Source.calib(Base, Target.PosX, Target.PosZ);
            Output += format("#{0}: {1:.4f} → ERR: {2:.4f} ± {3:.4f}\n", Source.Data.size(), 180/pi * Error, 180/pi * Source.Emean, 180/pi * Source.Esigma);
            Output += format("/tp {0:.2f} 240.00 {1:.2f}\n", Target.PosX + 60 * cos(Angle), Target.PosZ + 60 * sin(Angle));
        }
        else if (not Source.Data.empty()) {
            const Stronghold Target(Base, Source);
            for (const auto& Point : Target.Data)
                Output += format("{0:>6.3f}% → ({1:.0f}, {2:.0f})\n", 100 * Point.Prob, Point.PosX, Point.PosZ);
            Output += format("({0:.0f} ± {1:.0f}, {2:.0f} ± {3:.0f})\n", Target.Xmean, Target.Xsigma, Target.Zmean, Target.Zsigma);
        }
        return Output;
    }

};