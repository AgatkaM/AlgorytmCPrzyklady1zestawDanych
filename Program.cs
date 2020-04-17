using System;

using static funkcje.MetodaThomasa;

namespace ConsoleApp1
{
    class Program
    {
        static void Main(string[] args)
        {
            double cS=1000, cL=1290, roS=2679, roL=2380, L=390000, tempKrzep=930, u0=1013,uI=200,lambdaL=104,lambdaS=240;
            int lbIteracji = 10;
            double a,stalaA,stalaSigma;
            double entalpiaL = 3000000000, entalpiaS = 250000000;
            double b = 10,tGwiazdka=10000;
            double alfa=0.5;
            
            double[] temp0= { 200,200,200,200,200}; 
            int n = temp0.Length;
            double[] aVec = new double[n];
            double[] bVec = new double[n];
            double[] cVec = new double[n];
            double[] dVec = new double[n];
            double[,] tabH=new double[lbIteracji+1, n];
            double[,] tabHKreska = new double[lbIteracji+1, n];
            double[,] tabHDaszek = new double[lbIteracji+1, n];
            double[,] tabH2Kreska = new double[lbIteracji+1, n];
            double[,] tabH2Daszek = new double[lbIteracji + 1, n];
            int[] tabTx = new int[n]; /*tablica, w której zapisane są czasy tx topnienia dla danego punktu x*/ 
            double deltaX = b / n,deltaT=tGwiazdka/lbIteracji;
            double suma;
           

            funkcje.MetodaThomasa thomas = new funkcje.MetodaThomasa(); /* obiekt klasy MetodaThomasa do rozwiązania układu równań */
            double[] wynikThomas;

            funkcje.Temperatura tempZEntalpii = new funkcje.Temperatura(); /* obiekt klasy Temperatura do wyznaczenia temperatury na podstawie entalpii */
            double temp;
            
            funkcje.FunkcjaG funkcjaG = new funkcje.FunkcjaG(); /* obiekt klasy funkcjaG */
            double wynikG;


            for (int i = 0; i < n; i++)
                tabTx[i] = lbIteracji; /* ustawiamy czas topnienia dla punktu x na lbIteracji */

            for (int i = 0; i < temp0.Length; i++)
            {
                if (temp0[i] > tempKrzep)
                {
                    tabH[0, i] = cS * roS * tempKrzep + cL * roL * (temp0[i] - tempKrzep) + L * roS; /* Etap 1: wyznaczamy entalpie w chwili 0 */
                }
                else
                    tabH[0, i] = cS * roS * temp0[i];
                }

            /*Console.Write("tabH0: ");
            for (int i = 0; i < n; i++)
                Console.Write("{0} ", tabH[0, i]);*/

            stalaSigma = 1 / (SpecialFunction.gamma(0.5) * (1 - alfa) * Math.Pow(deltaT, alfa)); /* wartość funkcji sigma */

            for (int j = 1; j < lbIteracji; j++)
            {
                a = lambdaL / (cL*roL); /* Etap 2a: sprowadzamy całość do fazy ciekłej */
                stalaA = -a / (Math.Pow(deltaX, 2));
                

                for (int i = 0; i < temp0.Length; i++)
                    tabHKreska[j-1, i] = Math.Max(tabH[j-1, i], entalpiaL);


                aVec[0] = 0;
                for (int i = 1; i < n - 1; i++)
                    aVec[i] = stalaA;
                aVec[n-1] = 2*stalaA;

                bVec[0] = 1;
                for (int i = 1; i < n; i++)
                    bVec[i] = stalaSigma-2*stalaA;

                cVec[0] = 0;
                for (int i = 1; i < n - 1; i++)
                    cVec[i] = stalaA;
                cVec[n-1] = 0;

                dVec[0] = cS * roS * (tempKrzep) + cL * roL * (u0 - tempKrzep) + L * roS; /* obliczenie entalpii dla brzegowej temperatury u0 */



                for (int i = 1; i < n; i++)
                {
                    suma = 0;
                    for (int it = 2; it <= j; it++)
                    {
                        suma = suma + (Math.Pow(it, 1 - alfa) - Math.Pow(it - 1, 1 - alfa)) * (tabHKreska[j-it+1, i]-tabHKreska[j-it,i]);
                    }

                    dVec[i] = stalaSigma * tabHKreska[j-1, i] - stalaSigma * suma+ cL*roL*funkcjaG.g(i, tabTx[i], tabHKreska, j, stalaSigma, alfa, cL, cS, roL, roS);
                }

                wynikThomas= thomas.nowa(aVec, bVec, cVec, dVec); /* rozwiązanie układu równań metodą Thomasa */

                for (int i = 0; i < n; i++)
                    tabHDaszek[j, i] = wynikThomas[i];

                /*Console.Write("tabHDaszek: ");
                for (int i = 0; i < n; i++)
                    Console.Write("{0} ", tabHDaszek[j, i]);*/


                for (int i = 0; i < n; i++)
                    tabHKreska[j,i]=tabHDaszek[j, i]+(tabH[j-1,i]-tabHKreska[j-1,i]); /* korekta entalpii */

                /*Console.Write("tabHKreska: ");
                for (int i = 0; i < n; i++)
                    Console.Write("{0} ", tabHKreska[j, i]);*/


                a = lambdaS / (cS * roS); /* Etap 2b: sprowadzamy całość do fazy stałej  */

                for (int i = 0; i < n; i++)
                    tabH2Kreska[j-1, i] = Math.Min(tabHDaszek[j, i], entalpiaS);

                aVec[0] = 0;
                for (int i = 1; i < n - 1; i++)
                    aVec[i] = stalaA;
                aVec[n - 1] = 0;

                
                for (int i = 0; i < n-1; i++)
                    bVec[i] = stalaSigma - 2 * stalaA;
                bVec[n-1] = 1;

                cVec[0] = 2*stalaA;
                for (int i = 1; i < n - 1; i++)
                    cVec[i] = stalaA;
                cVec[n - 1] = 0;


                for (int i = 0; i < n - 1; i++)
                {
                    suma = 0;
                    for (int it = 2; it <= j; it++)
                    {
                        suma = suma + (Math.Pow(it, 1 - alfa) - Math.Pow(it - 1, 1 - alfa)) * (tabH2Kreska[j - it + 1, i] - tabH2Kreska[j - it, i]);
                    }

                    dVec[i] = stalaSigma * tabH2Kreska[j-1, i] - stalaSigma * suma + cS * roS * funkcjaG.g(i, tabTx[i], tabH2Kreska, j, stalaSigma, alfa, cL, cS, roL, roS); ;
                }

                dVec[n - 1] = cS * roS * uI;

                wynikThomas = thomas.nowa(aVec, bVec, cVec, dVec); /* rozwiązanie układu równań metodą Thomasa */


                for (int i = 0; i < n; i++)
                    tabH2Daszek[j , i] = wynikThomas[i];

                /*Console.Write("tabH2Daszek: ");
                for (int i = 0; i < n; i++)
                    Console.Write("{0} ", tabH2Daszek[j, i]);*/

                for (int i = 0; i < n; i++)
                    tabH[j, i] = tabH2Daszek[j, i] + (tabHKreska[j, i] - tabH2Kreska[j-1, i]);

                Console.Write("tabH: ");
                for (int i = 0; i < n; i++)
                {
                    Console.Write("{0} ", tabH[j, i]);
                }

                Console.WriteLine(" ");

                for(int i=0;i<n;i++)
                {
                    if(tempZEntalpii.temperatura(tabH[j,i])==tempKrzep&&tempZEntalpii.temperatura(tabH[j-1,i])!=tempKrzep)
                    {
                    tabTx[i]=j;
                    }

                }
            }

            Console.Write("Wartość tabTx:");
            for(int i=0;i<n;i++)
            {
                
                Console.Write("{0} ", tabTx[i]);
            }

            Console.WriteLine("");

            temp = tempZEntalpii.temperatura(3000000000);
            Console.WriteLine("Wartosc temperatury: {0}",temp);

            wynikG = funkcjaG.g(n-1, tabTx[n-1], tabH, 3, stalaSigma, alfa, cL, cS, roL, roS);
            Console.WriteLine("Wartość g: {0} ",wynikG);

            Console.WriteLine("Wartość funkcji gamma: {0} ", SpecialFunction.gamma(0.5));


            Console.ReadKey();
            return;
        }
    }
}
