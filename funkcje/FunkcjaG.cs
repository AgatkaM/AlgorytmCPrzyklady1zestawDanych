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
            double tempKrzep = 1773,entalpiaL=10693590000,entalpiaS=8563590000,L=270000;

            if (t > tx)
            {
                for (int it = 1; it <= tx; it++)
                {
                    suma = suma + (Math.Pow(it, 1 - alfa) - Math.Pow(it - 1, 1 - alfa)) * (tempZEntalpii.temperatura(tab[tx-it+1,x], tempKrzep, cL, roL, cS, roS, L, entalpiaL, entalpiaS) -tempZEntalpii.temperatura(tab[tx-it,x], tempKrzep, cL, roL, cS, roS, L, entalpiaL, entalpiaS));
                }

                g = ((cL * roL) / (cS * roS) - 1) * stalaSig * suma;
            }
            else
                g = 0;

            return g;
        }
    }

}
  
