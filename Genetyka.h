#pragma once
#include <string>
#include <vector>

using namespace std;

const int LICZBA_KODONOW = 64;
const int LICZBA_AMINOKWASOW = 20;

struct WynikORF
{
    int numer;
    int pozycjaStart;
    int pozycjaStop;
    vector<string> aminokwasy;
};

vector<string> WczytajPlik(string sciezka);
string PolaczLinieDNA(vector<string> linie);
string DNAnaMRNA(string dna);

bool CzyStart(string kodon);
bool CzyStop(string kodon);

string KodonNaAminokwas(string kodon);

vector<WynikORF> AnalizujMRNA(string mrna);

void WypiszWyniki(string dna, string mrna, vector<WynikORF> wyniki);