export module iTrace;
export import std;
export import "cubiomes/finders.h";
export import "cubiomes/util.h";
using namespace std;
using namespace numbers;

export class iTrace {

protected:

	struct Constants {
		struct ring { double Rmin, Rmax, Count, *Distr; };
		double Chunk, PosGen, PosMid; vector<ring> data;
		static double table7[], table12[], table16[];
		Constants(long long Base) {
			if (Base < MC_1_9) {
				Chunk = 16, PosGen = 4, PosMid = 4;
				data.emplace_back(460, 1332, 3, table7);
			}
			else if (Base < MC_1_13) {
				Chunk = 16, PosGen = 4, PosMid = 0;
				data = {
					{  1228,  2868,  3, table12 },
					{  4300,  5940,  6, table12 },
					{  7372,  9012, 10, table12 },
					{ 10444, 12084, 15, table12 },
					{ 13516, 15156, 21, table12 },
					{ 16588, 18228, 28, table12 },
					{ 19660, 21300, 36, table12 },
					{ 22732, 24372,  9, table12 }
				};
			}
			else if (Base < MC_1_19) {
				Chunk = 16, PosGen = 4, PosMid = 0;
				data = {
					{  1228,  2868,  3, table16 },
					{  4300,  5940,  6, table16 },
					{  7372,  9012, 10, table16 },
					{ 10444, 12084, 15, table16 },
					{ 13516, 15156, 21, table16 },
					{ 16588, 18228, 28, table16 },
					{ 19660, 21300, 36, table16 },
					{ 22732, 24372,  9, table16 }
				};
			}
			else {
				Chunk = 16, PosGen = 4, PosMid = 8;
				data = {
					{  1228,  2868,  3, table16 },
					{  4300,  5940,  6, table16 },
					{  7372,  9012, 10, table16 },
					{ 10444, 12084, 15, table16 },
					{ 13516, 15156, 21, table16 },
					{ 16588, 18228, 28, table16 },
					{ 19660, 21300, 36, table16 },
					{ 22732, 24372,  9, table16 }
				};
			}
		}
	};

	struct Endereyes {
		struct eye { double PosX, PosZ, Yaw, Range; };
		mutable double Emean, Esigma; vector<eye> data;
		double calib(const Constants& Base, double PosX, double PosZ) const {
			PosX = Base.Chunk * (floor(PosX / Base.Chunk) + 0.5) - Base.PosMid;
			PosZ = Base.Chunk * (floor(PosZ / Base.Chunk) + 0.5) - Base.PosMid;
			double Error = 0, Esum1 = 0, Esum2 = 0, Rsum2 = 0;
			for (const auto& eye : data) {
				Error = remainder(eye.Yaw - atan2(PosZ - eye.PosZ, PosX - eye.PosX), 2 * pi);
				Esum1 += Error / data.size(), Esum2 += Error * Error / (data.size() - 1);
				Rsum2 += eye.Range * eye.Range / data.size();
			}
			Esum2 = fdim(Esum2, Esum1 * Esum1 * data.size() / (data.size() - 1));
			Emean = Esum1, Esigma = sqrt(fdim(Esum2, Rsum2 / 3));
			return Error;
		}
		double solve(const Constants& Base, double PosX, double PosZ) const {
			PosX = Base.Chunk * (floor(PosX / Base.Chunk) + 0.5) - Base.PosMid;
			PosZ = Base.Chunk * (floor(PosZ / Base.Chunk) + 0.5) - Base.PosMid;
			double Prob = 0, Radius = hypot(PosX + Base.PosMid, PosZ + Base.PosMid);
			for (const auto& ring : Base.data)
				if (Radius > ring.Rmin and Radius < ring.Rmax)
					Prob = ring.Count * (ring.Rmax - ring.Rmin) / (Base.Chunk * Radius);
			if (Prob == 0) return Prob;
			double Dmin = +numeric_limits<double>::infinity();
			double Dmax = -numeric_limits<double>::infinity();
			for (const auto& eye : data) {
				Dmin = fmin(Dmin, hypot(eye.PosX + Base.PosMid, eye.PosZ + Base.PosMid) - hypot(PosX - eye.PosX, PosZ - eye.PosZ));
				Dmax = fmax(Dmax, hypot(eye.PosX + Base.PosMid, eye.PosZ + Base.PosMid) + hypot(PosX - eye.PosX, PosZ - eye.PosZ));
				double Error = remainder(eye.Yaw - atan2(PosZ - eye.PosZ, PosX - eye.PosX), 2 * pi);
				double Emin = (Error - Emean - eye.Range) / sqrt2;
				double Emax = (Error - Emean + eye.Range) / sqrt2;
				Prob *= Esigma > 0 ? Esigma * (erf(Emax / Esigma) - erf(Emin / Esigma)) / eye.Range : Emin < 0 and Emax > 0;
			}
			vector<pair<double, double>> cache;
			for (const auto& ring : Base.data) {
				auto Distr = [&ring](double Radius) { return Radius < ring.Rmax ? Radius < ring.Rmin ? 0 : ring.Distr[size_t(Radius - ring.Rmin)] : 1; };
				if (Radius > ring.Rmin and Radius < ring.Rmax) {
					Prob *= Distr(Radius + 0.5 * Base.Chunk) - Distr(Radius - 0.5 * Base.Chunk);
					for (double Index = 1; Index < ring.Count; Index++) {
						double Angle = 2 * pi * Index / ring.Count + atan2(PosZ, PosX);
						for (const auto& eye : data) {
							double Scale = (PosX + Base.PosMid) * (eye.PosX + Base.PosMid) + (PosZ + Base.PosMid) * (eye.PosZ + Base.PosMid);
							double Coeff = (eye.PosX + Base.PosMid) * cos(Angle) + (eye.PosZ + Base.PosMid) * sin(Angle);
							double Delta = Radius * Radius + Coeff * Coeff - 2 * Scale;
							if (Delta > 0) cache.emplace_back(Coeff - sqrt(Delta), Coeff + sqrt(Delta));
						}
						ranges::sort(cache, ranges::less());
						double Pmult = 1, Rmax = -numeric_limits<double>::infinity();
						for (const auto& line : cache) {
							Pmult -= Distr(fmax(Rmax, line.second)) - Distr(fmax(Rmax, line.first));
							Rmax = fmax(Rmax, line.second);
						}
						Prob *= Pmult, cache.clear();
					}
				}
				else if (Dmax > ring.Rmin and Dmin < ring.Rmax) {
					double Psum = 0, Count = round(pi * (ring.Rmin + ring.Rmax) / (Base.Chunk * ring.Count));
					for (double Step = 0.5; Step < Count; Step++) {
						double Prob = 1 / Count;
						for (double Index = 0; Index < ring.Count; Index++) {
							double Angle = 2 * pi * (Index + Step / Count) / ring.Count;
							for (const auto& eye : data) {
								double Scale = (PosX + Base.PosMid) * (eye.PosX + Base.PosMid) + (PosZ + Base.PosMid) * (eye.PosZ + Base.PosMid);
								double Coeff = (eye.PosX + Base.PosMid) * cos(Angle) + (eye.PosZ + Base.PosMid) * sin(Angle);
								double Delta = Radius * Radius + Coeff * Coeff - 2 * Scale;
								if (Delta > 0) cache.emplace_back(Coeff - sqrt(Delta), Coeff + sqrt(Delta));
							}
							ranges::sort(cache, ranges::less());
							double Pmult = 1, Rmax = -numeric_limits<double>::infinity();
							for (const auto& pair : cache) {
								Pmult -= Distr(fmax(Rmax, pair.second)) - Distr(fmax(Rmax, pair.first));
								Rmax = fmax(Rmax, pair.second);
							}
							Prob *= Pmult, cache.clear();
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
		struct str { double PosX, PosZ, Prob; };
		double Xmean, Xsigma, Zmean, Zsigma; vector<str> data;
		Stronghold(long long Base, long long Seed) {
			thread_local Generator Source;
			thread_local StrongholdIter Target;
			setupGenerator(&Source, Base, false);
			applySeed(&Source, DIM_OVERWORLD, Seed);
			initFirstStronghold(&Target, Base, Seed);
			nextStronghold(&Target, &Source);
			Xmean = Target.pos.x, Xsigma = 0;
			Zmean = Target.pos.z, Zsigma = 0;
			data.emplace_back(Target.pos.x, Target.pos.z, 1);
		}
		Stronghold(const Constants& Base, const Endereyes& Source) {
			vector<size_t> step; vector<vector<pair<double, double>>> cache;
			for (const auto& eye : Source.data) {
				step.emplace_back(), cache.emplace_back();
				double Radius = hypot(eye.PosX + Base.PosMid, eye.PosZ + Base.PosMid);
				double Dmin = +numeric_limits<double>::infinity();
				double Dmax = +numeric_limits<double>::infinity();
				for (const auto& ring : Base.data) {
					double Angle = pi / ring.Count;
					Dmin = fmin(Dmin, Radius < ring.Rmax ? Radius < ring.Rmin ? ring.Rmin - Radius : 0 : Radius - ring.Rmax);
					Dmax = fmin(Dmax, hypot(Radius * sin(Angle), fmax(Radius * cos(Angle) - ring.Rmin, ring.Rmax - Radius * cos(Angle))));
				}
				double Angle = eye.Yaw - Source.Emean;
				double Amin = Angle - eye.Range - 4 * Source.Esigma;
				double Amax = Angle + eye.Range + 4 * Source.Esigma;
				double PosBox[] = { eye.PosX + Dmin * cos(Amin), eye.PosX + Dmax * cos(Amin), eye.PosX + Dmin * cos(Amax), eye.PosX + Dmax * cos(Amax) };
				double Xmin = Base.Chunk * (round((ranges::min(PosBox) + Base.PosMid) / Base.Chunk) + 0.5) - Base.PosMid;
				double Xmax = Base.Chunk * (round((ranges::max(PosBox) + Base.PosMid) / Base.Chunk) + 0.5) - Base.PosMid;
				for (double PosX = Xmin; PosX < Xmax; PosX += Base.Chunk) {
					double Zmin = +numeric_limits<double>::infinity();
					double Zmax = -numeric_limits<double>::infinity();
					if (PosX > fmin(PosBox[0], PosBox[1]) and PosX < fmax(PosBox[0], PosBox[1])) {
						Zmin = fmin(Zmin, eye.PosZ + (PosX - eye.PosX) * tan(Amin));
						Zmax = fmax(Zmax, eye.PosZ + (PosX - eye.PosX) * tan(Amin));
					}
					if (PosX > fmin(PosBox[2], PosBox[3]) and PosX < fmax(PosBox[2], PosBox[3])) {
						Zmin = fmin(Zmin, eye.PosZ + (PosX - eye.PosX) * tan(Amax));
						Zmax = fmax(Zmax, eye.PosZ + (PosX - eye.PosX) * tan(Amax));
					}
					if (PosX > fmin(PosBox[0], PosBox[2]) and PosX < fmax(PosBox[0], PosBox[2])) {
						Zmin = fmin(Zmin, eye.PosZ + Dmin / sin(Angle) - (PosX - eye.PosX) / tan(Angle));
						Zmax = fmax(Zmax, eye.PosZ + Dmin / sin(Angle) - (PosX - eye.PosX) / tan(Angle));
					}
					if (PosX > fmin(PosBox[1], PosBox[3]) and PosX < fmax(PosBox[1], PosBox[3])) {
						Zmin = fmin(Zmin, eye.PosZ + Dmax / sin(Angle) - (PosX - eye.PosX) / tan(Angle));
						Zmax = fmax(Zmax, eye.PosZ + Dmax / sin(Angle) - (PosX - eye.PosX) / tan(Angle));
					}
					Zmin = Base.Chunk * (round((Zmin + Base.PosMid) / Base.Chunk) + 0.5) - Base.PosMid;
					Zmax = Base.Chunk * (round((Zmax + Base.PosMid) / Base.Chunk) + 0.5) - Base.PosMid;
					for (double PosZ = Zmin; PosZ < Zmax; PosZ += Base.Chunk)
						cache.back().emplace_back(PosX, PosZ);
				}
			}
			double Xsum1 = 0, Xsum2 = 0, Zsum1 = 0, Zsum2 = 0;
			if (not cache.empty()) {
				step.emplace_back(), cache.emplace_back();
				double Psum = 0; LOOP: size_t IDmin = 0, IDmax = 0;
				for (size_t Index = 0; Index + 1 < cache.size(); Index++) {
					if (step[Index] == cache[Index].size()) goto EXIT;
					else if (cache[Index][step[Index]] < cache[IDmin][step[IDmin]]) step[Index]++, IDmin = Index;
					else if (cache[Index][step[Index]] < cache[IDmax][step[IDmax]]) step[Index]++;
					else if (cache[Index][step[Index]] > cache[IDmax][step[IDmax]]) step[IDmax]++, IDmax = Index;
				}
				if (IDmin == IDmax) {
					cache.back().emplace_back(cache[0][step[0]]);
					Psum += Source.solve(Base, cache[0][step[0]].first, cache[0][step[0]].second);
					for (size_t Index = 0; Index + 1 < cache.size(); step[Index++]++);
				}
				goto LOOP; EXIT:
				for (const auto& pair : cache.back()) {
					double PosX = Base.Chunk * floor(pair.first / Base.Chunk) + Base.PosGen;
					double PosZ = Base.Chunk * floor(pair.second / Base.Chunk) + Base.PosGen;
					double Prob = Source.solve(Base, pair.first, pair.second) / Psum;
					Xsum1 += PosX * Prob, Xsum2 += PosX * PosX * Prob;
					Zsum1 += PosZ * Prob, Zsum2 += PosZ * PosZ * Prob;
					if (Prob > 0) data.emplace_back(PosX, PosZ, Prob);
				}
			}
			Xmean = Xsum1, Xsigma = sqrt(fdim(Xsum2, Xsum1 * Xsum1));
			Zmean = Zsum1, Zsigma = sqrt(fdim(Zsum2, Zsum1 * Zsum1));
			ranges::sort(data, ranges::greater(), &str::Prob);
		}
	};

private:

	long long Base = MC_NEWEST, Seed = 0;

	Endereyes Source = { pi/720, pi/3600 };

public:

	string operator()(const string& Input) {
		static regex Pattern[] = {
			regex("CHECK", regex::icase),
			regex("CLEAR", regex::icase),
			regex("VER (\\S+)", regex::icase),
			regex("CAL (\\S+)", regex::icase),
			regex("ERR (\\S+) (\\S+)", regex::icase),
			regex("ADD (\\S+) (\\S+) (\\S+)", regex::icase),
			regex("/execute in (?:minecraft:)?overworld run tp @s (\\S+) (\\S+) (\\S+) (\\S+) (\\S+)", regex::icase)
		};
		thread_local default_random_engine RNG;
		size_t Index; smatch Value; string Output;
		for (Index = 0; Index < size(Pattern); Index++)
			if (regex_match(Input, Value, Pattern[Index])) break;
		switch (Index) {
		case 0: Index = 0; break;
		case 1: Source.data.clear(); break;
		case 2: Base = str2mc(Value[1].str().data()); break;
		case 3: Seed = stoll(Value[1]), RNG.seed(stoll(Value[1])); break;
		case 4: Source.Emean = pi/180 * stod(Value[1]), Source.Esigma = pi/180 * stod(Value[2]); break;
		case 5: Source.data.emplace_back(stod(Value[1]), stod(Value[2]), pi/180 * (stod(Value[3]) + 90), pi/3600); break;
		case 6: Source.data.emplace_back(stod(Value[1]), stod(Value[3]), pi/180 * (stod(Value[4]) + 90), pi/36000); break;
		default: throw invalid_argument("invalid iTrace input");
		}
		if (Index == 0) {
			Output += format("VER: {0}  SEED: {1}  ERR: {2:.4f} ± {3:.4f}\n", mc2str(Base), Seed, 180/pi * Source.Emean, 180/pi * Source.Esigma);
			for (const auto& eye : Source.data)
				Output += format("#{0}: ({1:.1f}, {2:.1f}, {3:.2f} ± {4:.1g})\n", ++Index, eye.PosX, eye.PosZ, 180/pi * eye.Yaw - 90, 180/pi * eye.Range);
		}
		else if (Seed != 0) {
			const auto& str = Stronghold(Base, Seed).data.back();
			double Angle = uniform_real_distribution(-pi, pi)(RNG);
			double Error = Source.calib(Base, str.PosX, str.PosZ);
			Output += format("#{0}: {1:.4f} → ERR: {2:.4f} ± {3:.4f}\n", Source.data.size(), 180/pi * Error, 180/pi * Source.Emean, 180/pi * Source.Esigma);
			Output += format("/tp {0:.2f} 240.00 {1:.2f}\n", str.PosX + 60 * cos(Angle), str.PosZ + 60 * sin(Angle));
		}
		else if (not Source.data.empty()) {
			Stronghold Target(Base, Source);
			for (const auto& str : Target.data)
				Output += format("{0:>6.3f}% → ({1:.0f}, {2:.0f})\n", 100 * str.Prob, str.PosX, str.PosZ);
			Output += format("({0:.0f} ± {1:.0f}, {2:.0f} ± {3:.0f})\n", Target.Xmean, Target.Xsigma, Target.Zmean, Target.Zsigma);
		}
		return Output;
	}

};