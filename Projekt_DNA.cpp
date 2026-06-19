#include <iostream>
#include <string>
#include <vector>
#include "Genetyka.h"

using namespace std;

int main()
{
    string sciezka;

    cout << "Podaj sciezke do pliku txt z DNA." << endl;
    cout << "Jesli plik nazywa sie dna.txt i jest w folderze projektu, nacisnij ENTER: ";

    getline(cin, sciezka);

    if (sciezka == "")
    {
        sciezka = "dna.txt";
    }

    if (sciezka.length() >= 2 && sciezka[0] == '"' && sciezka[sciezka.length() - 1] == '"')
    {
        sciezka = sciezka.substr(1, sciezka.length() - 2);
    }

    vector<string> linie = WczytajPlik(sciezka);

    if (linie.size() == 0)
    {
        return 1;
    }

    string dna = PolaczLinieDNA(linie);
    string mrna = DNAnaMRNA(dna);

    vector<WynikORF> wyniki = AnalizujMRNA(mrna);

    WypiszWyniki(dna, mrna, wyniki);

    return 0;
}