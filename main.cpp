import iTrace;
using namespace std;

int main() {
	iTrace Instance; string Input;
	ofstream save("history.log", ios::app);
	for (ifstream load("history.log"); getline(load, Input).good(); Instance(Input));
	cout << Instance("CHECK") << endl;
	while (getline(cin, Input).good() and not Input.empty())
		cout << Instance(Input) << endl, save << Input << endl;
}