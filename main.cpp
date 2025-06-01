import iTrace;
using namespace std;

int main() try {
	iTrace Instance; string Input;
	for (ifstream file("history.log"); getline(file, Input); Instance(Input));
	cout << Instance("CHECK") << endl;
	ofstream("README.MD") << iTrace::Intro << endl;
	while (getline(cin, Input) and not Input.empty()) {
		cout << Instance(Input) << endl;
		ofstream("history.log", ios::app) << Input << endl;
	}
}
catch (const exception& err) {
	cerr << err.what() << endl;
	return main();
}