using System;

namespace funkcje
{

    class Temperatura
    {
        double temp;
        public double temperatura(double entalpia,double temperaturaKrzepniecia, double cL,double roL,double cS, double roS, double L, double entalpiaL, double entalpiaS)
        {

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
  /*  class Program
    {
        static void Main(string[] args)
        {
            double temp;

            Temperatura t = new Temperatura();
            temp = t.temperatura(2500000000);

            
            Console.WriteLine(temp);
        }
    }*/
}
