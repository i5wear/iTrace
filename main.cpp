import iTrace;
using namespace std;
using namespace chrono;

int main() try {
	iTrace Instance; string Input;
	for (ifstream file("history.log"); getline(file, Input); Instance(Input));
	cout << Instance("CHECK") << endl;
	ofstream("README.MD") << iTrace::Intro << endl;
	while (getline(cin, Input) and not Input.empty()) {
		auto begin = high_resolution_clock::now();
		string Output = Instance(Input);
		auto end = high_resolution_clock::now();
		Output += format("{0}\n", duration_cast<microseconds>(end - begin));
		cout << Output << endl;
		ofstream("history.log", ios::app) << Input << endl;
	}
}
catch (const exception& err) {
	cerr << err.what() << endl;
	return main();
}