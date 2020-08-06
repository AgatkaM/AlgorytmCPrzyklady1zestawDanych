using System;

namespace funkcje
{

    class Entalpia
    {
        double temp;
        public double entalpia(double temperatura,double temperaturaKrzepniecia, double cL,double roL,double cS, double roS, double L)
        {
            double ent;

            if (temperatura > temperaturaKrzepniecia)
            {
                ent = cL * roL * (temperatura - temperaturaKrzepniecia) + cS * roS * temperaturaKrzepniecia + L * roS;
            } 
            else 
            {
                ent = cS * roS * temperatura;
            }
           

            return ent;
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
