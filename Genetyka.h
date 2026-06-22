#pragma once
#include <string>

using namespace std;

const int LICZBA_KODONOW = 64;
const int LICZBA_AMINOKWASOW = 20;
const int MAKS_ORF = 5000;

struct WynikORF
{
    int numer;
    int pozycjaStart;
    int pozycjaStop;
    int liczbaAminokwasow;
};

string WczytajPlik(string sciezka);
string DNAnaMRNA(string dna);

bool CzyStart(string kodon);
bool CzyStop(string kodon);

string KodonNaAminokwas(string kodon);

int AnalizujMRNA(string mrna, WynikORF wyniki[]);

void WypiszWyniki(string dna, string mrna, WynikORF wyniki[], int liczbaWynikow);