#include "Genetyka.h"
#include <iostream>
#include <fstream>

using namespace std;

const string KODONY[LICZBA_KODONOW] = {
    "UUU", "UUC", "UUA", "UUG",
    "UCU", "UCC", "UCA", "UCG",
    "UAU", "UAC", "UAA", "UAG",
    "UGU", "UGC", "UGA", "UGG",

    "CUU", "CUC", "CUA", "CUG",
    "CCU", "CCC", "CCA", "CCG",
    "CAU", "CAC", "CAA", "CAG",
    "CGU", "CGC", "CGA", "CGG",

    "AUU", "AUC", "AUA", "AUG",
    "ACU", "ACC", "ACA", "ACG",
    "AAU", "AAC", "AAA", "AAG",
    "AGU", "AGC", "AGA", "AGG",

    "GUU", "GUC", "GUA", "GUG",
    "GCU", "GCC", "GCA", "GCG",
    "GAU", "GAC", "GAA", "GAG",
    "GGU", "GGC", "GGA", "GGG"
};

const string AMINOKWASY[LICZBA_KODONOW] = {
    "Fenyloalanina", "Fenyloalanina", "Leucyna", "Leucyna",
    "Seryna", "Seryna", "Seryna", "Seryna",
    "Tyrozyna", "Tyrozyna", "STOP", "STOP",
    "Cysteina", "Cysteina", "STOP", "Tryptofan",

    "Leucyna", "Leucyna", "Leucyna", "Leucyna",
    "Prolina", "Prolina", "Prolina", "Prolina",
    "Histydyna", "Histydyna", "Glutamina", "Glutamina",
    "Arginina", "Arginina", "Arginina", "Arginina",

    "Izoleucyna", "Izoleucyna", "Izoleucyna", "Metionina",
    "Treonina", "Treonina", "Treonina", "Treonina",
    "Asparagina", "Asparagina", "Lizyna", "Lizyna",
    "Seryna", "Seryna", "Arginina", "Arginina",

    "Walina", "Walina", "Walina", "Walina",
    "Alanina", "Alanina", "Alanina", "Alanina",
    "Kwas asparaginowy", "Kwas asparaginowy", "Kwas glutaminowy", "Kwas glutaminowy",
    "Glicyna", "Glicyna", "Glicyna", "Glicyna"
};

const string LISTA_AMINOKWASOW[LICZBA_AMINOKWASOW] = {
    "Fenyloalanina", "Leucyna", "Seryna", "Tyrozyna", "Cysteina",
    "Tryptofan", "Prolina", "Histydyna", "Glutamina", "Arginina",
    "Izoleucyna", "Metionina", "Treonina", "Asparagina", "Lizyna",
    "Walina", "Alanina", "Kwas asparaginowy", "Kwas glutaminowy", "Glicyna"
};

vector<string> WczytajPlik(string sciezka)
{
    vector<string> linie;
    ifstream plik(sciezka);

    if (!plik.is_open())
    {
        cerr << "Nie mozna otworzyc pliku: " << sciezka << endl;
        return linie;
    }

    string linia;

    while (getline(plik, linia))
    {
        linie.push_back(linia);
    }

    plik.close();

    return linie;
}

string PolaczLinieDNA(vector<string> linie)
{
    string dna = "";

    for (int i = 0; i < linie.size(); i++)
    {
        for (int j = 0; j < linie[i].length(); j++)
        {
            char znak = linie[i][j];

            if (znak == 'a') znak = 'A';
            if (znak == 't') znak = 'T';
            if (znak == 'g') znak = 'G';
            if (znak == 'c') znak = 'C';

            if (znak == 'A' || znak == 'T' || znak == 'G' || znak == 'C')
            {
                dna += znak;
            }
        }
    }

    return dna;
}

string DNAnaMRNA(string dna)
{
    string mrna = "";

    for (int i = 0; i < dna.length(); i++)
    {
        if (dna[i] == 'T')
        {
            mrna += 'U';
        }
        else
        {
            mrna += dna[i];
        }
    }

    return mrna;
}

bool CzyStart(string kodon)
{
    return kodon == "AUG";
}

bool CzyStop(string kodon)
{
    return kodon == "UAA" || kodon == "UAG" || kodon == "UGA";
}

string KodonNaAminokwas(string kodon)
{
    for (int i = 0; i < LICZBA_KODONOW; i++)
    {
        if (KODONY[i] == kodon)
        {
            return AMINOKWASY[i];
        }
    }

    return "NIEZNANY";
}

vector<WynikORF> AnalizujMRNA(string mrna)
{
    vector<WynikORF> wyniki;

    int numer = 1;
    int i = 0;
    int dlugosc = mrna.length();

    while (i <= dlugosc - 3)
    {
        string kodon = mrna.substr(i, 3);

        if (CzyStart(kodon))
        {
            WynikORF wynik;

            wynik.numer = numer;
            wynik.pozycjaStart = i;
            wynik.pozycjaStop = -1;

            bool znalezionoStop = false;

            for (int j = i; j <= dlugosc - 3; j = j + 3)
            {
                string kodonAktualny = mrna.substr(j, 3);

                if (CzyStop(kodonAktualny))
                {
                    wynik.pozycjaStop = j;
                    znalezionoStop = true;
                    break;
                }

                string aminokwas = KodonNaAminokwas(kodonAktualny);

                if (aminokwas != "NIEZNANY")
                {
                    wynik.aminokwasy.push_back(aminokwas);
                }
            }

            if (znalezionoStop)
            {
                wyniki.push_back(wynik);
                numer++;

                i = wynik.pozycjaStop + 3;
            }
            else
            {
                i++;
            }
        }
        else
        {
            i++;
        }
    }

    return wyniki;
}

void WypiszWyniki(string dna, string mrna, vector<WynikORF> wyniki)
{
    cout << "===== WYNIK ANALIZY DNA =====" << endl;
    cout << "Dlugosc DNA: " << dna.length() << endl;
    cout << "DNA:  " << dna << endl;
    cout << "mRNA: " << mrna << endl;
    cout << endl;

    if (wyniki.size() == 0)
    {
        cout << "Nie znaleziono kompletnej sekwencji START-STOP." << endl;
        return;
    }

    int ilosci[LICZBA_AMINOKWASOW] = { 0 };

    cout << "Znalezione sekwencje START-STOP: " << wyniki.size() << endl;
    cout << "Pozycje sa liczone od 0 w sekwencji mRNA." << endl;
    cout << endl;

    for (int i = 0; i < wyniki.size(); i++)
    {
        cout << "ORF " << wyniki[i].numer << endl;
        cout << "START: " << wyniki[i].pozycjaStart << endl;
        cout << "STOP:  " << wyniki[i].pozycjaStop << endl;

        cout << "Aminokwasy: ";

        for (int j = 0; j < wyniki[i].aminokwasy.size(); j++)
        {
            cout << wyniki[i].aminokwasy[j];

            if (j < wyniki[i].aminokwasy.size() - 1)
            {
                cout << ", ";
            }

            for (int k = 0; k < LICZBA_AMINOKWASOW; k++)
            {
                if (wyniki[i].aminokwasy[j] == LISTA_AMINOKWASOW[k])
                {
                    ilosci[k]++;
                }
            }
        }

        cout << endl;
        cout << endl;
    }

    cout << "===== LISTA WYKRYTYCH AMINOKWASOW I ICH ILOSC =====" << endl;

    for (int i = 0; i < LICZBA_AMINOKWASOW; i++)
    {
        if (ilosci[i] > 0)
        {
            cout << LISTA_AMINOKWASOW[i] << ": " << ilosci[i] << endl;
        }
    }
}