using System;

namespace funkcje
{
    class MetodaThomasa
    {
        public double[] nowa(double[] tabA, double[] tabB, double[] tabC, double[] tabD)
        {
            int n = tabA.Length;
            double m;
            //float[] tabA = { 0, 3, 1, 1, 2 };
            //float[] tabB = { 2, 1, 2, 1, 2 };
            //float[] tabC = { 2, 1, 4, 1, 0 };
            //float[] tabD = { 2, 6, 4, 1, 7 };
            double[] tabBeta = new double[n];
            double[] tabGamma = new double[n];
            double[] tabX = new double[n];

            tabBeta[0] = -tabC[0] / tabB[0];
            tabGamma[0] = tabD[0] / tabB[0];

            for (int i = 1; i < n; i++)
            {
                m = tabA[i] * tabBeta[i - 1] + tabB[i];
                tabBeta[i] = -tabC[i] / m;
                tabGamma[i] = (tabD[i] - tabA[i] * tabGamma[i - 1]) / m;
            }

            tabX[n - 1] = tabGamma[n - 1];

            for (int i = n - 2; i >= 0; i--)
            {
                tabX[i] = tabBeta[i] * tabX[i + 1] + tabGamma[i];
            }

            /*for (int i = 0; i < n; i++)
                Console.Write("{0} ", tabX[i]);*/

            return tabX;
        }
    }

  /*  class Program
    {
        static void Main(string[] args)
        {
            double[] tabX0;
            double[] tabA = { 0, 3, 1, 1, 2 };
            double[] tabB = { 2, 1, 2, 1, 2 };
            double[] tabC = { 2, 1, 4, 1, 0 };
            double[] tabD = { 2, 6, 4, 1, 7 };


            MetodaThomasa t = new MetodaThomasa();
          
            tabX0 = t.nowa(tabA,tabB,tabC,tabD);
            Console.Write(tabX0[0]);

            Console.ReadKey();

        }
    }*/
}
