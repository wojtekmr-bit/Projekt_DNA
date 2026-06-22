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

string WczytajPlik(string sciezka)
{
    string dna = "";
    ifstream plik(sciezka);

    if (!plik.is_open())
    {
        cerr << "Nie mozna otworzyc pliku: " << sciezka << endl;
        return dna;
    }

    string linia;

    while (getline(plik, linia))
    {
        for (int i = 0; i < linia.length(); i++)
        {
            char znak = linia[i];

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

    plik.close();

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

int AnalizujMRNA(string mrna, WynikORF wyniki[])
{
    int liczbaWynikow = 0;
    int numer = 1;
    int i = 0;
    int dlugosc = mrna.length();

    while (i <= dlugosc - 3 && liczbaWynikow < MAKS_ORF)
    {
        string kodon = mrna.substr(i, 3);

        if (CzyStart(kodon))
        {
            int pozycjaStop = -1;

            for (int j = i; j <= dlugosc - 3; j = j + 3)
            {
                string kodonAktualny = mrna.substr(j, 3);

                if (CzyStop(kodonAktualny))
                {
                    pozycjaStop = j;
                    break;
                }
            }

            if (pozycjaStop != -1)
            {
                wyniki[liczbaWynikow].numer = numer;
                wyniki[liczbaWynikow].pozycjaStart = i;
                wyniki[liczbaWynikow].pozycjaStop = pozycjaStop;
                wyniki[liczbaWynikow].liczbaAminokwasow = (pozycjaStop - i) / 3;

                liczbaWynikow++;
                numer++;

                i = pozycjaStop + 3;
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

    return liczbaWynikow;
}

void WypiszWyniki(string dna, string mrna, WynikORF wyniki[], int liczbaWynikow)
{
    cout << "===== WYNIK ANALIZY DNA =====" << endl;
    cout << "Dlugosc DNA: " << dna.length() << endl;
    cout << "DNA:  " << dna << endl;
    cout << "mRNA: " << mrna << endl;
    cout << endl;

    if (liczbaWynikow == 0)
    {
        cout << "Nie znaleziono kompletnej sekwencji START-STOP." << endl;
        return;
    }

    int ilosci[LICZBA_AMINOKWASOW] = { 0 };

    cout << "Znalezione sekwencje START-STOP: " << liczbaWynikow << endl;
    cout << "Pozycje sa liczone od 0 w sekwencji mRNA." << endl;
    cout << endl;

    for (int i = 0; i < liczbaWynikow; i++)
    {
        cout << "ORF " << wyniki[i].numer << endl;
        cout << "START: " << wyniki[i].pozycjaStart << endl;
        cout << "STOP:  " << wyniki[i].pozycjaStop << endl;
        cout << "Dlugosc sekwencji aminokwasow: " << wyniki[i].liczbaAminokwasow << endl;

        cout << "Aminokwasy: ";

        for (int j = wyniki[i].pozycjaStart; j < wyniki[i].pozycjaStop; j = j + 3)
        {
            string kodon = mrna.substr(j, 3);
            string aminokwas = KodonNaAminokwas(kodon);

            cout << aminokwas;

            if (j + 3 < wyniki[i].pozycjaStop)
            {
                cout << ", ";
            }

            for (int k = 0; k < LICZBA_AMINOKWASOW; k++)
            {
                if (aminokwas == LISTA_AMINOKWASOW[k])
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