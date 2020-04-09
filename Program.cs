using System;

namespace ConsoleApp1
{
    class Program
    {
        static void Main(string[] args)
        {
            double cS=1000, cL=1290, roS=2679, roL=2380, L=390000, tempKrzep=930;
            int k = 1;
            double entalpiaL = 3000000000, entalpiaS = 250000000;
            double[] temp0= { 1013,900,500};
            double[,] tabH=new double[k,temp0.Length];
            double[,] tabHKreska = new double[k,temp0.Length];

            for (int i = 0; i < temp0.Length; i++)
                tabH[0,i] = cS*roS*tempKrzep+cL*roL*(temp0[i]-tempKrzep)+L*roS;

            for (int i = 0; i < temp0.Length; i++)
                tabHKreska[0,i] = Math.Max(tabH[0,i], entalpiaL);

            for (int i=0; i<temp0.Length; i++)
                Console.Write("{0} ", tabHKreska[0,i]);
            Console.WriteLine("");


            for (int i = 0; i < temp0.Length; i++)
                tabHKreska[0,i] = Math.Max(tabH[0,i], entalpiaL);


            Console.ReadKey();
            return;
        }
    }
}
