#include <iostream>
#include <string>
#include "Genetyka.h"

using namespace std;

int main()
{
    string dna = WczytajPlik("dna.txt");

    if (dna.length() == 0)
    {
        return 1;
    }

    string mrna = DNAnaMRNA(dna);

    WynikORF wyniki[MAKS_ORF];

    int liczbaWynikow = AnalizujMRNA(mrna, wyniki);

    WypiszWyniki(dna, mrna, wyniki, liczbaWynikow);

    return 0;
}