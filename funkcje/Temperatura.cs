using System;

namespace ConsoleApp3
{

    class Temperatura
    {
        double temp;
        public double temperatura(double entalpia)
        {
            double temperaturaKrzepniecia=930, cL=1290, roL=2380, cS=1000, roS=2679, L=390000,entalpiaL=3000000000,entalpiaS=250000000;

            if (entalpia > entalpiaL)
            {
                temp = temperaturaKrzepniecia + 1 / (cL * roL) * (entalpia - cS * roS * temperaturaKrzepniecia - L * roS);
            } 
            else if (entalpia < entalpiaS)
            {
                temp = 1 / (cS * roS) * entalpia;
            }
            else
                temp = temperaturaKrzepniecia;
      

            return temp;
        }

    }
    class Program
    {
        static void Main(string[] args)
        {
            double temp;

            Temperatura t = new Temperatura();
            temp = t.temperatura(2500000000);

            
            Console.WriteLine(temp);
        }
    }
}
