using System;

using static funkcje.MetodaThomasa;

namespace ConsoleApp1
{
    class Program
    {
        static void Main(string[] args)
        {
            double cS= 1000, cL= 1290, roS= 2679, roL= 2380, L= 390000, tempKrzep= 930, u0= 1213, uI= 298, lambdaL=104,lambdaS= 240;
            //double cS = 690, cL = 800, roS = 7500, roL = 7000, L = 270000, tempKrzep = 1773, u0 = 2213, uI = 323, lambdaL = 33, lambdaS = 30;
            //double cS = 10, cL = 10, roS = 25, roL = 25, L = 10, tempKrzep = 1, u0,uI, /*u0 = 1573, uI =973,*/ lambdaL = 250, lambdaS = 250;
            int lbIteracji = 1000;
            double a,stalaA,stalaSigma;
            //double entalpiaL = 3536280000, entalpiaS = 2491470000;
            //double entalpiaL = 11200275000, entalpiaS = 9175275000;
            //double entalpiaL = 500, entalpiaS = 250;
            double entalpiaS = cS*roS*tempKrzep, entalpiaL = entalpiaS+L*roS;
            double b = 2,tGwiazdka=100;
            double alfa=0.5;
            int granica = 0;
            int n = 100;
            double deltaX = b /n, deltaT = 1;


            //double[] temp0= { 50,50,50,50,50}; 
            //int n = temp0.Length;

           
            double[] temp0 = new double[n+1];
         

            for (int i=0;i<=n;i++)
            {
                //temp0[i] = tempKrzep-i*dT;
                temp0[i] = uI;
                //temp0[i] =Math.Exp(-deltaX*i); /*warunek początkowy T0 */
            }
            
            double[] aVec = new double[n+1];
            double[] bVec = new double[n+1];
            double[] cVec = new double[n+1];
            double[] dVec = new double[n+1];
            double[,] tabH=new double[lbIteracji+1, n+1];
            double[,] tabHKreska = new double[lbIteracji+1, n+1];
            double[,] tabHDaszek = new double[lbIteracji+1, n+1];
            double[,] tabH2Kreska = new double[lbIteracji+1, n+1];
            double[,] tabH2Daszek = new double[lbIteracji + 1, n+1];
            int[] tabTx = new int[n+1]; /*tablica, w której zapisane są czasy tx topnienia dla danego punktu x*/ 
            
            double suma=0;
            int j = 1;
            int l = 10;
            double[] temperatury = new double[n+1];


            funkcje.MetodaThomasa thomas = new funkcje.MetodaThomasa(); /* obiekt klasy MetodaThomasa do rozwiązania układu równań */
            double[] wynikThomas;

            funkcje.Temperatura tempZEntalpii = new funkcje.Temperatura(); /* obiekt klasy Temperatura do wyznaczenia temperatury na podstawie entalpii */
            double temp;
            
            funkcje.FunkcjaG funkcjaG = new funkcje.FunkcjaG(); /* obiekt klasy funkcjaG */
            double wynikG;

            funkcje.Entalpia entZtemp = new funkcje.Entalpia();
            double ent;

            


            for (int i = 0; i <= n; i++)
                tabTx[i] = lbIteracji; /* ustawiamy czas topnienia dla punktu x na lbIteracji */

            

            for (int i = 0; i < temp0.Length; i++)
            {
                if (temp0[i] > tempKrzep)
                {
                    tabH[0, i] = entZtemp.entalpia(temp0[i], tempKrzep, cL, roL, cS, roS, L); /* Etap 1: wyznaczamy entalpie w chwili 0 na podstawie T0*/
                }
                else
                    tabH[0, i] = entZtemp.entalpia(temp0[i], tempKrzep, cL, roL, cS, roS, L);
            }

            //Console.Write("tabH0[n-1]: {0}", tabH[0, 0]);

            /*Console.Write("tabH0: ");
            for (int i = 0; i < n; i++)
                Console.Write("{0} ", tabH[0, i]);*/
    

            stalaSigma = 1 / (SpecialFunction.gamma(1-alfa) * (1 - alfa) * Math.Pow(deltaT, alfa)); /* wartość funkcji sigma */



            //for (int j = 1; j < lbIteracji; j++)
            while (granica < l&&j<=lbIteracji)
            {
                a = lambdaL / (cL*roL); /* Etap 2a: sprowadzamy całość do fazy ciekłej */
                //a = 1;
                stalaA = a / (Math.Pow(deltaX, 2));
                

                for (int i = 0; i < temp0.Length; i++)
                    tabHKreska[j-1, i] = Math.Max(tabH[j-1, i], entalpiaL);

                //Console.WriteLine(" ");
                //Console.Write("tabHKreska: ");
                //for (int i = 0; i < n; i++)
                //    Console.Write("{0} ", tabHKreska[j-1, i]);

                //Console.Write("tabHKreska[n-1]: {0}", tabHKreska[j-1, 0]);

                aVec[0] = 0;
                for (int i = 1; i < n; i++)
                    aVec[i] = -stalaA;
                aVec[n] = -2*stalaA;

                bVec[0] = 1;
                for (int i = 1; i < n+1; i++)
                    bVec[i] = stalaSigma+2*stalaA;

                cVec[0] = 0;
                for (int i = 1; i < n ; i++)
                    cVec[i] = -stalaA;
                cVec[n] = 0;

                //u0 = Math.Exp(a*(j)*deltaT);
                dVec[0] = entZtemp.entalpia(u0, tempKrzep, cL, roL, cS, roS, L); /* obliczenie entalpii dla brzegowej temperatury u0 */
                Console.WriteLine(dVec[0]>entalpiaL);

                for (int i = 1; i <=n; i++)
                {
                    suma = 0;
                    
                    for (int it = 2; it <= j; it++)
                    {
                        suma = suma + (Math.Pow(it, 1 - alfa) - Math.Pow(it - 1, 1 - alfa)) * (tabHKreska[j-it+1, i]-tabHKreska[j-it,i]);
                    }

                    dVec[i] = stalaSigma * tabHKreska[j-1, i] - stalaSigma * suma+ cL*roL*funkcjaG.g(i, tabTx[i], tabHKreska, j, stalaSigma, alfa, cL, cS, roL, roS, tempKrzep, entalpiaL, entalpiaS, L);
                }
                

                wynikThomas= thomas.nowa(aVec, bVec, cVec, dVec); /* rozwiązanie układu równań metodą Thomasa */

                for (int i = 0; i <=n; i++)
                    tabHDaszek[j, i] = wynikThomas[i];

               //Console.Write("tabHDaszek[n-1]: {0}", tabHDaszek[j, 0]);

                //Console.WriteLine(" ");
                //Console.Write("tabHDaszek wyznaczone z Thomasa: ");
                //for (int i = 0; i < n; i++)
                //    Console.Write("{0} ", tabHDaszek[j, i]);


                for (int i = 0; i <=n; i++)
                    tabHKreska[j,i]=tabHDaszek[j, i]+(tabH[j-1,i]-tabHKreska[j-1,i]); /* korekta entalpii */

                //Console.Write("tabHKreska[n-1]: {0}", tabHKreska[j, 0]);

                /* Console.WriteLine(" ");
                 Console.Write("tabHKreska po korekcie: ");
                 for (int i = 0; i < n; i++)
                     Console.Write("{0} ", tabHKreska[j, i]); */


                a = lambdaS / (cS * roS); /* Etap 2b: sprowadzamy całość do fazy stałej  */
                //a = 1;
                stalaA = a / (Math.Pow(deltaX, 2));


                for (int i = 0; i <= n; i++)
                {
                    tabH2Kreska[j - 1, i] = Math.Min(tabHKreska[j, i], entalpiaS);
                }

                //Console.Write("tabH2Kreska[n-1]: {0}", tabH2Kreska[j-1, 0]);

                //Console.WriteLine(" ");
                //Console.Write("tabH2Kreska (minimum): ");
                //for (int i = 0; i < n; i++)
                //    Console.Write("{0} ", tabH2Kreska[j-1, i]);

                aVec[0] = 0;
                for (int i = 1; i < n ; i++)
                    aVec[i] = -stalaA;
                aVec[n] = 0;

                
                for (int i = 0; i < n; i++)
                    bVec[i] = stalaSigma + 2 * stalaA;
                bVec[n] = 1;

                cVec[0] = -2*stalaA;
                for (int i = 1; i < n; i++)
                    cVec[i] = -stalaA;
                cVec[n] = 0;


                for (int i = 0; i <n; i++)
                {
                    suma = 0;
                    for (int it = 2; it <= j; it++)
                    {
                        suma = suma + (Math.Pow(it, 1 - alfa) - Math.Pow(it - 1, 1 - alfa)) * (tabH2Kreska[j - it + 1, i] - tabH2Kreska[j - it, i]);
                    }
                    
                    dVec[i] = stalaSigma * tabH2Kreska[j-1, i] - stalaSigma * suma + cL * roL * funkcjaG.g(i, tabTx[i], tabH2Kreska, j, stalaSigma, alfa, cL, cS, roL, roS, tempKrzep, entalpiaL, entalpiaS, L); ;
                }
                //uI = Math.Exp(a*(j)*deltaT-2);
                dVec[n] = entZtemp.entalpia(uI, tempKrzep, cL, roL, cS, roS, L);
                Console.WriteLine(uI);

                wynikThomas = thomas.nowa(aVec, bVec, cVec, dVec); /* rozwiązanie układu równań metodą Thomasa */


                for (int i = 0; i <=n; i++)
                    tabH2Daszek[j , i] = wynikThomas[i];

                //Console.Write("tabH2Daszek[n-1]: {0}", tabH2Daszek[j, 0]);

                //Console.WriteLine(" ");
                // Console.Write("tabH2Daszek z thomasa: ");
                // for (int i = 0; i < n; i++)
                // Console.Write("{0} ", tabH2Daszek[j, i]);

                for (int i = 0; i <=n; i++)
                    tabH[j, i] = tabH2Daszek[j, i] + (tabHKreska[j, i] - tabH2Kreska[j-1, i]); /* korekta */

                //Console.Write("tabH[n-1]: {0}", tabH[j, 0]);

                //Console.WriteLine(" ");
                //Console.Write("tabH po korekcie: ");
                //for (int i = 0; i < n; i++)
                //{
                //    Console.Write("{0} ", tabH[j, i]);
                //}

                //Console.WriteLine(" ");
                //Console.Write("temp tabH: ");
                //for (int i = 0; i < n; i++)
                //{
                //    Console.Write("{0} ", tempZEntalpii.temperatura(tabH[j, i]));
                //}

                //Console.WriteLine(" ");

                for (int i=0;i<=n;i++)
                {
                    if(tabH[j,i]>entalpiaS&&tabH[j-1,i]<=entalpiaS)
                    {
                    tabTx[i]=j;
                    }

                }

                for (int i=1;i<n;i++)
                {
                    if(tabH[j,i]<=entalpiaS&&tabH[j,i+1]<=entalpiaS&&tabH[j,i-1]>entalpiaS)
                    {
                        granica = i;
                    }
                }


                Console.WriteLine("{0}", j);
                j++;
            }

            Console.Write("tabH: ");
            for (int i = 0; i <=n; i++)
            {
                Console.Write("{0} ", tabH[j-1, i]);
            }

            string path = @"C:\Users\agata.chmielowska\Desktop\AlgorytmC-\entalpie.csv";
            System.IO.StreamWriter str = new System.IO.StreamWriter(path, false);

            for (int i = 0; i <=n; i++)
            {
                str.WriteLine("{0} ", tabH[j-1,i]);
            }
            str.Close();

            Console.WriteLine("");
            Console.Write("temp tabH: ");
            for (int i = 0; i <=n; i++)
            {
                temperatury[i] = tempZEntalpii.temperatura(tabH[j - 1, i], tempKrzep, cL, roL, cS, roS, L, entalpiaL, entalpiaS);
                Console.Write("{0} ", temperatury[i]);
            }

            path = @"C:\Users\agata.chmielowska\Desktop\AlgorytmC-\temperatury.csv";
            str = new System.IO.StreamWriter(path, false);

            for (int i = 0; i <=n; i++)
            {
                str.WriteLine("{0} ", temperatury[i]);
            }
            str.Close();

            Console.WriteLine("");
            Console.Write("Wartość tabTx:");
            for(int i=0;i<=n;i++)
            {
                
                Console.Write("{0} ", tabTx[i]);
            }

            path = @"C:\Users\agata.chmielowska\Desktop\AlgorytmC-\granica.csv";
            str = new System.IO.StreamWriter(path, false);

            for (int i = 0; i <=n; i++)
            {
                str.WriteLine("{0} ", tabTx[i]*deltaT);
            }
            str.Close();

            Console.WriteLine("");
            Console.WriteLine("Program zakończył działanie po {0} iteracjach", j-1);

            Console.WriteLine("");
            Console.WriteLine("Granica znajduje sie w punkcie: {0}", granica);

            temp = tempZEntalpii.temperatura(entalpiaL, tempKrzep, cL, roL, cS, roS, L, entalpiaL, entalpiaS);
            Console.WriteLine("Wartosc temperatury: {0}",temp);

            wynikG = funkcjaG.g(5, tabTx[5], tabH, 20, stalaSigma, alfa, cL, cS, roL, roS, tempKrzep, entalpiaL, entalpiaS, L);
            Console.WriteLine("Wartość g: {0} ",wynikG);

            Console.WriteLine("Wartość funkcji gamma: {0} ", SpecialFunction.gamma(0.5));





            Console.ReadKey();
            return;
        }
    }
}
