import iTrace;
using namespace std;

int main() {
	iTrace Instance; string Input;
	for (ifstream load("history.log"); getline(load, Input).good(); Instance(Input));
	cout << Instance("CHECK") << endl;
	ofstream save("history.log", ios::app);
	while (getline(cin, Input).good() and not Input.empty())
		cout << Instance(Input) << endl, save << Input << endl;
}