Für die drei Reaktionsmodelle existieren drei Skripte:
	Chemical_Degradation_final.ipynb
	Electrochemical_Degradation_final.ipynb
	Cuppled_Degradation_final.ipynb

Alle drei Skripte sind identisch aufgebaut:
	Bibliotheken bzw. Funktionen werden importiert (z. B. rate_const_H2O2)
	Berechnungen werden durchgeführt, die sich nicht mehr ändern (Anzahl reaktiver Zentren, Fluorid-Konzentration etc.)
	Definition von Konstanten und Funktionen (rxn, gas_trans, deg_rate, H2_price)
	Ausführung der Simulation unter verschiedenen Parametern in unterschiedlichen Darstellungsformen

Für die Ausführung der Simulation wird die Funktion deg_rate ausgeführt. Diese nimmt Temperatur T, Stromdichte curden, Kathodendruck P_cat, Anodendruck P_an und Membrandicke h_memb als Eingangsparameter. Je nach Skript nimmt sie auch andere Parameter an, die für den jeweiligen Teil relevant sind (z. B. f_Fe, womit die Eisenflussrate im elektrochemischen Teil eingestellt werden kann, was für den chemischen Teil nicht von Bedeutung ist). Dabei ist zu beachten, dass deg_rate immer nur einzelne Werte annimmt, nie Listen oder Arrays! Soll über einer Liste eines Parameters simuliert werden, so muss dies mit einem Loop gelöst werden.
Nach Eingabe der Parameter berechnet deg_rate alle davon abhängigen Größen und führt das Reaktionsnetzwerk rxn aus. Ausgabe von deg_rate ist eine Liste verschiedener Größen, an deren nullter Stelle i. d. R. die attack_rate in der Einheit mol/(l*s) steht.
Cuppled_degradation_final.ipynb enthält auch die Funktion H2_price. Diese nimmt Temperatur T, Stromdichte curden, Kathodendruck P_cat, Anodendruck P_an, Membrandicke h_memb, aktive Membranfläche A, Lebenszeit lifetime_h (in Stunden !) und Energiepreis energy_price (in €/Wh !) als Eingangsparameter. Ausgegeben wird eine Liste mit dem H2-Preis pro Kilogramm, der Gesamtmasse an produziertem H2, CAPEX, OPEX und Energiekosten.


Weitere Skripte sind:
	H2O2_formation.ipynb; Dieses simuliert und validiert das Zellmodell. Am Ende steht die Funktion rate_const_H2O2, welche Temperatur T, Kathodendruck P_cat, Anodendruck P_an und Membrandicke h_memb als Eingangsparameter verwendet. Ausgabe ist eine Liste, an deren nullter Stelle eine weitere Liste steht. Diese innere Liste enthält die Geschwindigkeitskonstante der H2O2-Bildung in Abhängigkeit von 10.000 Werten der Stromdichte von 0 bis 3 A/cm². Um die Geschwindigkeitskonstante bei einer Stromdichte zu erhalten, wird jener Wert der 10.000 ermittelt, der dem Eingangswert am nächsten liegt. Der Index dieses Näherungswertes wird verwendet, um die dort liegende Geschwindigkeitskonstante zu ermitteln. > k=rate_const_H2O2(T,P_cat,P_an,h_memb)[0][abs(j-curden).argmin()]
	gas_trans_v2.ipynb; Dieses simuliert den Massentransport durch die Membran nach Schalenbach. Es enthält die Funktion gas_trans, die in den Modellen verwendet wird.


An verschiedenen Stellen werden Daten aus Veröffentlichungen zum Vergleich herangezogen. Diese Daten sind in folgenden Textdateien hinterlegt:
	Chandesris_FER.txt; Diese enthält x-Werte [0:6] und y-Werte [6:12] der FER vs. curden bei 60 °C sowie x-Werte [12:15] und y-Werte [15:18] der FER vs. curden bei 80 °C. Abb. 3 aus DOI: 10.1016/j.ijhydene.2014.11.111
	Hellsing_desorption.txt; Diese enthält x-Werte [0:19] und y-Werte [19:38] der OH-Desorptionsraten sowie x-Werte [38:67] und y-Werte [67:96] der H2O-Desorptionsraten. Abb. 1 aus DOI: 10.1016/0021-9517(91)90258-6
	Laconti_FER.txt; Diese enthält x-Werte [0:6] und y-Werte [6:12] der FER vs. T sowie zwei x-Werte [12:14] und zwei y-Werte [14:16] für die Ausgleichsgerade. Abb. 9 aus DOI: 10.1149/1.2214554
	Schalenbach_crossover.txt; Für H2 in O2 gegen Stromdichte enthält diese x-Werte [0:10] und y-Werte [10:20] bei Gleichdruck und x-Werte [20:29] und y-Werte [29:38] bei Differentialdruck. Abb. 3 aus DOI: 10.1016/j.ijhydene.2013.09.013