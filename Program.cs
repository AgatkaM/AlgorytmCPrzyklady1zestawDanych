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
            double b = 10,tGwiazdka=100;
            double alfa=0.5;
            
            double[] temp0= { 1013,900,500};
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
            double[] tabTx = new double[n];
            double deltaX = b / n,deltaT=tGwiazdka/lbIteracji;
            double suma = 0;

            for (int i = 0; i < n; i++)
                tabTx[i] = tGwiazdka;

            for (int i = 0; i < temp0.Length; i++)
                tabH[0,i] = cS*roS*tempKrzep+cL*roL*(temp0[i]-tempKrzep)+L*roS;

            for (int j = 0; j < lbIteracji; j++)
            {
                a = lambdaL / (cL*roL);
                stalaA = -a / (Math.Pow(deltaX, 2));
                stalaSigma = 1 / 0.5 *(1-alfa)*Math.Pow(deltaT,alfa);

                for (int i = 0; i < temp0.Length; i++)
                    tabHKreska[j, i] = Math.Max(tabH[j, i], entalpiaL);


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

                dVec[0] = cS * roS * (tempKrzep) + cL * roL * (u0 - tempKrzep) + L * roS;



                for (int i = 1; i < n; i++)
                {
                    for (int it = 2; it <= j; it++)
                    {
                        suma = suma + (Math.Pow(it, 1 - alfa) - Math.Pow(it - 1, 1 - alfa)) * (tabHKreska[j-it+1, i]-tabHKreska[j-it,i]);
                    }

                    dVec[i] = stalaSigma * tabHKreska[j, i] - stalaSigma * suma;
                }

                //dVec

                for (int i = 0; i < n; i++)
                    tabHDaszek[j+1, i] = tabHKreska[j, i];

                for (int i = 0; i < n; i++)
                    tabHKreska[j+1,i]=tabHDaszek[j+1, i]+(tabH[j,i]-tabHKreska[j,i]);

                for (int i = 0; i < n; i++)
                    Console.Write("{0} ", tabHKreska[j, i]);


                a = lambdaS / (cS * roS);
                for (int i = 0; i < n; i++)
                    tabH2Kreska[j, i] = Math.Min(tabHDaszek[j+1, i], entalpiaS);

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
                    for (int it = 2; it <= j; it++)
                    {
                        suma = suma + (Math.Pow(it, 1 - alfa) - Math.Pow(it - 1, 1 - alfa)) * (tabH2Kreska[j - it + 1, i] - tabH2Kreska[j - it, i]);
                    }

                    dVec[i] = stalaSigma * tabHKreska[j, i] - stalaSigma * suma;
                }

                dVec[n - 1] = cS * roS * uI;

                //dVec

                for (int i = 0; i < n; i++)
                    tabH2Daszek[j + 1, i] = tabH2Kreska[j, i];

                for (int i = 0; i < n; i++)
                    tabH[j + 1, i] = tabH2Daszek[j + 1, i] + (tabHKreska[j+1, i] - tabH2Kreska[j, i]);

                for (int i = 0; i < n; i++)
                    Console.Write("{0} ", tabH[j+1, i]);

                Console.WriteLine("");

                for(int i=0;i<n;i++)
                {
                    //if(temperatura(tabH[j+1,i])==temperaturaKrzepniecia&&temperatura(tabH[j+1,i])!=temperaturaKrzepniecia)
                    //{
                    //tabTx[i]=j+1;
                    //}

                }
            }

            double tabX0;
            float[] tabA = { 0, 3, 1, 1, 2 };
            float[] tabB = { 2, 1, 2, 1, 2 };
            float[] tabC = { 2, 1, 4, 1, 0 };
            float[] tabD = { 2, 6, 4, 1, 7 };


            funkcje.MetodaThomasa t = new funkcje.MetodaThomasa();

            tabX0 = t.nowa(tabA, tabB, tabC, tabD);
            Console.Write(tabX0);

            Console.ReadKey();

            Console.ReadKey();
            return;
        }
    }
}
