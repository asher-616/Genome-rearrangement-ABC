#include <fstream>
#include <iostream>

int main() {
	std::fstream f;
	f.open("test.txt", std::fstream::in);
	std::cout << "managed to open";
	f.close();
}