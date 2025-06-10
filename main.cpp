import iTrace;
using namespace std;
using namespace chrono;

int main() {
	iTrace Instance; string Input;
	for (ifstream load{"history.log"}; getline(load, Input).good(); Instance(Input));
	cout << Instance("CHECK") << endl;
	ofstream("README.TXT", ios::noreplace) << iTrace::Intro << endl;
	while (getline(cin, Input).good() and not Input.empty()) {
		auto begin = high_resolution_clock::now();
		string Output = Instance(Input);
		auto end = high_resolution_clock::now();
		cout << Output << format("{0}\n", duration_cast<microseconds>(end - begin)) << endl;
		ofstream("history.log", ios::app) << Input << endl;
	}
}