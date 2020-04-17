using System;

namespace funkcje
{
    class FunkcjaG
    {
        public double g(int x, int tx, double[,] tab, int t, double stalaSig, double alfa, double cL, double cS, double roL, double roS)
        {
            funkcje.Temperatura tempZEntalpii = new funkcje.Temperatura();

            double g;
            double suma=0;

            if (t > tx)
            {
                for (int it = 1; it <= tx; it++)
                {
                    suma = suma + (Math.Pow(it, 1 - alfa) - Math.Pow(it - 1, 1 - alfa)) * (tempZEntalpii.temperatura(tab[tx-it+1,x])-tempZEntalpii.temperatura(tab[tx-it,x]));
                }

                g = ((cL * roL) / (cS * roS) - 1) * stalaSig * suma;
            }
            else
                g = 0;

            return g;
        }
    }

}
  
