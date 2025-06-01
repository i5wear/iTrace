export module iTrace;
export import std;
import "cubiomes/finders.h";
import "cubiomes/util.h";
using namespace std;
using namespace numbers;
#define VALUE " +([+-]?[0-9]+(?:[.][0-9]+)?)"

export class iTrace {

protected:
	
	struct Constants {
		struct ring { double Rmin, Rmax, Count, *Distr; };
		double chunk, PosGen, PosMid; vector<ring> data;
		static double table7[], table12[], table16[];
		Constants(MCVersion Base) {
			if (Base < MC_1_9) {
				chunk = 16, PosGen = 4, PosMid = 4;
				data = { {460, 1332, 3, table7} };
			}
			else if (Base < MC_1_13) {
				chunk = 16, PosGen = 4, PosMid = 0;
				data = {
					{ 1228,  2868,  3, table12},
					{ 4300,  5940,  6, table12},
					{ 7372,  9012, 10, table12},
					{10444, 12084, 15, table12},
					{13516, 15156, 21, table12},
					{16588, 18228, 28, table12},
					{19660, 21300, 36, table12},
					{22732, 24372,  9, table12}
				};
			}
			else if (Base < MC_1_19) {
				chunk = 16, PosGen = 4, PosMid = 0;
				data = {
					{ 1228,  2868,  3, table16},
					{ 4300,  5940,  6, table16},
					{ 7372,  9012, 10, table16},
					{10444, 12084, 15, table16},
					{13516, 15156, 21, table16},
					{16588, 18228, 28, table16},
					{19660, 21300, 36, table16},
					{22732, 24372,  9, table16}
				};
			}
			else {
				chunk = 16, PosGen = 4, PosMid = 8;
				data = {
					{ 1228,  2868,  3, table16},
					{ 4300,  5940,  6, table16},
					{ 7372,  9012, 10, table16},
					{10444, 12084, 15, table16},
					{13516, 15156, 21, table16},
					{16588, 18228, 28, table16},
					{19660, 21300, 36, table16},
					{22732, 24372,  9, table16}
				};
			}
		}
	};

	struct EnderEyes {
		struct eye { double PosX, PosZ, Yaw, Range; };
		double Emean, Evar; vector<eye> data;
		double calib(MCVersion Base, double PosX, double PosZ) {
			const Constants World(Base);
			PosX = World.chunk * (floor(PosX / World.chunk) + 0.5);
			PosZ = World.chunk * (floor(PosZ / World.chunk) + 0.5);
			double Error = 0; Emean = 0, Evar = 0;
			for (auto& eye : data) {
				eye.PosX += World.PosMid, eye.PosZ += World.PosMid;
				Error = remainder(eye.Yaw - atan2(PosZ - eye.PosZ, PosX - eye.PosX), 2 * pi);
				Emean += Error / data.size();
				Evar += Error * Error / (data.size() - 1);
				Evar -= eye.Range * eye.Range / (3 * data.size());
				eye.PosX -= World.PosMid, eye.PosZ -= World.PosMid;
			}
			Evar -= Emean * Emean;
			Evar = Evar > 0 ? sqrt(Evar) : 0;
			return Error;
		}
		double solve(MCVersion Base, double PosX, double PosZ) {
			const Constants World(Base); vector<pair<double, double>> cache;
			PosX = World.chunk * (floor(PosX / World.chunk) + 0.5);
			PosZ = World.chunk * (floor(PosZ / World.chunk) + 0.5);
			double Radius = hypot(PosX, PosZ), Prob = 0;
			for (const auto& ring : World.data)
				if (Radius > ring.Rmin and Radius < ring.Rmax)
					Prob = ring.Count * (ring.Rmax - ring.Rmin) / (Radius * World.chunk);
			if (not Prob) return Prob;
			double Dmin = +numeric_limits<double>::infinity();
			double Dmax = -numeric_limits<double>::infinity();
			for (auto& eye : data) {
				eye.PosX += World.PosMid, eye.PosZ += World.PosMid;
				Dmin = min(Dmin, hypot(eye.PosX, eye.PosZ) - hypot(PosX - eye.PosX, PosZ - eye.PosZ));
				Dmax = max(Dmax, hypot(eye.PosX, eye.PosZ) + hypot(PosX - eye.PosX, PosZ - eye.PosZ));
				if (Evar > 0) {
					double Error = remainder(eye.Yaw - atan2(PosZ - eye.PosZ, PosX - eye.PosX), 2 * pi);
					double Emin = (Error - Emean - eye.Range) / (sqrt2 * Evar);
					double Emax = (Error - Emean + eye.Range) / (sqrt2 * Evar);
					Prob *= Evar * (erf(Emax) - erf(Emin)) / eye.Range;
				}
			}
			for (const auto& ring : World.data) {
				double Count = round(pi * (ring.Rmin + ring.Rmax) / (World.chunk * ring.Count));
				auto Distr = [&ring](double Radius)
					{ return Radius < ring.Rmax ? Radius < ring.Rmin ? 0 : ring.Distr[size_t(Radius - ring.Rmin)] : 1; };
				if (Radius > ring.Rmin and Radius < ring.Rmax) {
					Prob *= Distr(Radius + 0.5 * World.chunk) - Distr(Radius - 0.5 * World.chunk);
					for (double Index = 1; Index < ring.Count; Index++) {
						double Angle = 2 * pi * Index / ring.Count + atan2(PosZ, PosX);
						for (const auto& eye : data) {
							double Coeff = eye.PosX * cos(Angle) + eye.PosZ * sin(Angle);
							double Delta = Radius * Radius + Coeff * Coeff - 2 * (PosX * eye.PosX + PosZ * eye.PosZ);
							if (Delta > 0) cache.emplace_back(Coeff - sqrt(Delta), Coeff + sqrt(Delta));
						}
						ranges::sort(cache);
						double Pmult = 1, Rmin = 0, Rmax = 0;
						for (const auto& line : cache) {
							Rmin = max(Rmax, line.first);
							Rmax = max(Rmax, line.second);
							Pmult -= Distr(Rmax) - Distr(Rmin);
						}
						Prob *= Pmult, cache.clear();
					}
				}
				else if (Dmax > ring.Rmin and Dmin < ring.Rmax) {
					double Psum = 0;
					for (double Step = 0.5; Step < Count; Step++) {
						double Prob = 1 / Count;
						for (double Index = 0; Index < ring.Count; Index++) {
							double Angle = 2 * pi * (Index + Step / Count) / ring.Count;
							for (const auto& eye : data) {
								double Coeff = eye.PosX * cos(Angle) + eye.PosZ * sin(Angle);
								double Delta = Radius * Radius + Coeff * Coeff - 2 * (PosX * eye.PosX + PosZ * eye.PosZ);
								if (Delta > 0) cache.emplace_back(Coeff - sqrt(Delta), Coeff + sqrt(Delta));
							}
							ranges::sort(cache);
							double Pmult = 1, Rmin = 0, Rmax = 0;
							for (const auto& line : cache) {
								Rmin = max(Rmax, line.first);
								Rmax = max(Rmax, line.second);
								Pmult -= Distr(Rmax) - Distr(Rmin);
							}
							Prob *= Pmult, cache.clear();
						}
						Psum += Prob;
					}
					Prob *= Psum;
				}
			}
			for (auto& eye : data)
				eye.PosX -= World.PosMid, eye.PosZ -= World.PosMid;
			return Prob;
		}
	};

	struct Stronghold {
		struct str { double PosX, PosZ, Prob; };
		double Xmean, Zmean, Xvar, Zvar; vector<str> data;
		Stronghold(MCVersion Base, long long Seed) {
			static Generator World;
			static StrongholdIter Target;
			setupGenerator(&World, Base, false);
			applySeed(&World, DIM_OVERWORLD, Seed);
			initFirstStronghold(&Target, Base, Seed);
			nextStronghold(&Target, &World);
			data.emplace_back(Target.pos.x, Target.pos.z, 1);
			Xmean = Target.pos.x, Zmean = Target.pos.z, Xvar = 0, Zvar = 0;
		}
		Stronghold(MCVersion Base, EnderEyes Source) {
			const Constants World(Base); vector<str> cache;
			auto order = [](const str& prev, const str& next)
				{ return prev.Prob != next.Prob ? prev.Prob > next.Prob : prev.PosX != next.PosX ? prev.PosX < next.PosX : prev.PosZ < next.PosZ; };
			for (auto& eye : Source.data) {
				eye.PosX += World.PosMid, eye.PosZ += World.PosMid;
				double Radius = hypot(eye.PosX, eye.PosZ);
				double Dmin = +numeric_limits<double>::infinity();
				double Dmax = +numeric_limits<double>::infinity();
				for (const auto& ring : World.data) {
					double Angle = pi / ring.Count;
					Dmin = min(Dmin, Radius < ring.Rmax ? Radius < ring.Rmin ? ring.Rmin - Radius : 0 : Radius - ring.Rmax);
					Dmax = min(Dmax, hypot(Radius * sin(Angle), max(Radius * cos(Angle) - ring.Rmin, ring.Rmax - Radius * cos(Angle))));
				}
				double Angle = eye.Yaw - Source.Emean;
				double Amin = Angle - eye.Range - 4 * Source.Evar;
				double Amax = Angle + eye.Range + 4 * Source.Evar;
				double PosBox[] = { eye.PosX + Dmin * cos(Amin), eye.PosX + Dmax * cos(Amin), eye.PosX + Dmin * cos(Amax), eye.PosX + Dmax * cos(Amax) };
				double Xmin = World.chunk * (round(ranges::min(PosBox) / World.chunk) + 0.5);
				double Xmax = World.chunk * (round(ranges::max(PosBox) / World.chunk) + 0.5);
				for (double PosX = Xmin; PosX < Xmax; PosX += World.chunk) {
					double Zmin = +numeric_limits<double>::infinity();
					double Zmax = -numeric_limits<double>::infinity();
					if (PosX > min(PosBox[0], PosBox[1]) and PosX < max(PosBox[0], PosBox[1])) {
						Zmin = min(Zmin, eye.PosZ + (PosX - eye.PosX) * tan(Amin));
						Zmax = max(Zmax, eye.PosZ + (PosX - eye.PosX) * tan(Amin));
					}
					if (PosX > min(PosBox[2], PosBox[3]) and PosX < max(PosBox[2], PosBox[3])) {
						Zmin = min(Zmin, eye.PosZ + (PosX - eye.PosX) * tan(Amax));
						Zmax = max(Zmax, eye.PosZ + (PosX - eye.PosX) * tan(Amax));
					}
					if (PosX > min(PosBox[0], PosBox[2]) and PosX < max(PosBox[0], PosBox[2])) {
						Zmin = min(Zmin, eye.PosZ + Dmin / sin(Angle) - (PosX - eye.PosX) / tan(Angle));
						Zmax = max(Zmax, eye.PosZ + Dmin / sin(Angle) - (PosX - eye.PosX) / tan(Angle));
					}
					if (PosX > min(PosBox[1], PosBox[3]) and PosX < max(PosBox[1], PosBox[3])) {
						Zmin = min(Zmin, eye.PosZ + Dmax / sin(Angle) - (PosX - eye.PosX) / tan(Angle));
						Zmax = max(Zmax, eye.PosZ + Dmax / sin(Angle) - (PosX - eye.PosX) / tan(Angle));
					}
					Zmin = World.chunk * (round(Zmin / World.chunk) + 0.5);
					Zmax = World.chunk * (round(Zmax / World.chunk) + 0.5);
					for (double PosZ = Zmin; PosZ < Zmax; PosZ += World.chunk)
						if (data.empty() or ranges::binary_search(data, str(PosX, PosZ, 0), order))
							cache.emplace_back(PosX, PosZ, 0);
				}
				data.swap(cache), cache.clear();
				eye.PosX -= World.PosMid, eye.PosZ -= World.PosMid;
				if (data.empty()) break;
			}
			double Psum = 0; Xmean = 0, Zmean = 0, Xvar = 0, Zvar = 0;
			for (auto& str : data) {
				str.PosX += World.PosGen - 0.5 * World.chunk;
				str.PosZ += World.PosGen - 0.5 * World.chunk;
				str.Prob = Source.solve(Base, str.PosX, str.PosZ);
				if (str.Prob > 0) {
					Psum += str.Prob;
					Xmean += str.Prob * str.PosX;
					Zmean += str.Prob * str.PosZ;
					cache.emplace_back(str.PosX, str.PosZ, str.Prob);
				}
			}
			Xmean /= Psum, Zmean /= Psum;
			data.swap(cache), ranges::sort(data, order);
			for (auto& str : data) {
				str.Prob /= Psum;
				Xvar += str.Prob * (str.PosX - Xmean) * (str.PosX - Xmean);
				Zvar += str.Prob * (str.PosZ - Zmean) * (str.PosZ - Zmean);
			}
			Xvar = sqrt(Xvar), Zvar = sqrt(Zvar);
		}
	};

private:

	MCVersion Base{ MC_1_16 };
	
	long long Seed{ false };

	EnderEyes Source{ 0, pi/18000 };

public:

	static constexpr char Intro[]{
		"# [iTrace](https://github.com/i5wear/iTrace)\n\n"
		"适用于多种场合的要塞计算器，由 i5wear 制作。\n\n"
		"基于准确的数学模型，其结果置信度不低于任何同类产品。\n\n"
		"本产品的研发和制作基于世界生成算法库 [Cubiomes](https://github.com/Cubitect/Cubiomes)。\n\n"
		"#### 使用帮助\n\n"
		"程序首次运行时，请先设定游戏版本，并设定误差或开始误差校准。\n\n"
		"对准末影之眼按 F3 + C，直接粘贴至本程序，即可快捷输入数据。\n\n"
		"若上述方式无法使用，则按下 F3，手动输入坐标和角度至本程序。\n\n"
		"正常创建一个世界，获取种子并粘贴至本程序，即可开始误差校准。\n\n"
		"终止本程序不会丢失任何数据，所有历史输入均被自动保存与读取。\n\n"
		"#### 命令列表\n\n"
		"添加末影之眼：ADD [X] [Z] [Yaw]\n\n"
		"设定游戏版本：VER 1.[N]\n\n"
		"设定测量误差：ERR [Mean] [SD]\n\n"
		"开始误差校准：CAL [Seed]\n\n"
		"终止误差校准：CAL 0\n\n"
		"查看所有数据：CHECK\n\n"
		"清除所有数据：CLEAR\n\n"
		"终止程序运行：空输入"
	};

	string operator()(const string& Input) {
		static regex Pattern[] = {
			regex("^ *CHECK *$", regex::icase),
			regex("^ *CLEAR *$", regex::icase),
			regex("^ *CAL" VALUE " *$", regex::icase),
			regex("^ *VER" VALUE " *$", regex::icase),
			regex("^ *ERR" VALUE VALUE " *$", regex::icase),
			regex("^ *ADD" VALUE VALUE VALUE " *$", regex::icase),
			regex("^ */execute in minecraft:overworld run tp @s" VALUE VALUE VALUE VALUE VALUE " *$", regex::icase)
		};
		size_t Index; smatch Value; string Output;
		for (Index = 0; Index < size(Pattern); Index++)
			if (regex_match(Input, Value, Pattern[Index])) break;
		switch (Index) {
		case 0: Index = 0; break;
		case 1: Source.data.clear(); break;
		case 2: Seed = stoll(Value[1]); break;
		case 3: Base = MCVersion(str2mc(Value[1].str().data())); break;
		case 4: Source.Emean = pi/180 * stod(Value[1]), Source.Evar = pi/180 * stod(Value[2]); break;
		case 5: Source.data.emplace_back(stod(Value[1]), stod(Value[2]), pi/180 * (stod(Value[3]) + 90), pi/3600); break;
		case 6: Source.data.emplace_back(stod(Value[1]), stod(Value[3]), pi/180 * (stod(Value[4]) + 90), pi/36000); break;
		default: return Output;
		}
		if (not Index) {
			Output += format("VER: {0}  Mean: {1:.4f}  SD: {2:.4f}\n", mc2str(Base), 180/pi * Source.Emean, 180/pi * Source.Evar);
			for (const auto& eye : Source.data)
				Output += format("#{0}: ({1:.1f}, {2:.1f}, {3:.2f} ± {4:.1g})\n", ++Index, eye.PosX, eye.PosZ, 180/pi * eye.Yaw - 90, 180/pi * eye.Range);
		}
		else if (Seed) {
			static default_random_engine RNG(time(nullptr));
			const auto& Target = Stronghold(Base, Seed).data.front();
			double Error = Source.calib(Base, Target.PosX, Target.PosZ);
			double Angle = uniform_real_distribution(-pi, pi)(RNG);
			if (not Source.data.empty())
				Output += format("ERR: {0:.4f}  Mean: {1:.4f}  SD: {2:.4f}\n", 180/pi * Error, 180/pi * Source.Emean, 180/pi * Source.Evar);
			Output += format("#{0}: /tp {1:.2f} 256 {2:.2f}\n", Source.data.size() + 1, Target.PosX + 64 * cos(Angle), Target.PosZ + 64 * sin(Angle));
		}
		else if (not Source.data.empty()) {
			Stronghold Target(Base, Source);
			for (const auto& str : Target.data)
				Output += format("{0:>6.3f}% -> ({1:.0f}, {2:.0f})\n", 100 * str.Prob, str.PosX, str.PosZ);
			if (not Target.data.empty())
				Output += format("({0:.0f} ± {1:.0f}, {2:.0f} ± {3:.0f})\n", Target.Xmean, Target.Xvar, Target.Zmean, Target.Zvar);
		}
		return Output;
	}

};