using System;
using System.Diagnostics;
using System.Threading;


namespace Com_Methods
{
    class CONST
    {
        //точность решения
        public static double EPS = 1e-15;
    }

    class Program
    {
        static void Main()
        {
            try
            {
                //прямые методы: метод Гаусса, LU-разложение, QR-разложение

                int N = 200;
                var A = new Matrix(N, N);
                var X_true = new Vector(N);

                //заполнение СЛАУ
                for (int i = 0; i < N; i++) 
                {
                    for (int j = 0; j < N; j++)
                    {
                        if (i != j) A.Elem[i][j] = 2 + 0.01 * i - 0.2 * j;
                        else A.Elem[i][j] = 100;
                    }
                    X_true.Elem[i] = 1;
                }

                //правая часть
                var F = A * X_true;

                //решатель
                var Solver = new Gauss_Method();
                //var Solver = new LU_Decomposition(A);
                //var Solver = new QR_Decomposition(A, QR_Decomposition.QR_Algorithm.Householder);

                //старт счетчика времени
                Stopwatch stopWatch = new Stopwatch();
                stopWatch.Start();
                
                //решение СЛАУ
                var X = Solver.Start_Solver(A,F); // для гаусс метода
                //var X = Solver.Start_Solver(F); //для всех остальных методов
                
                //стоп счетчика
                stopWatch.Stop();
                TimeSpan ts = stopWatch.Elapsed;

                //Вывод времени в миллисекундах
                string elapsedTime = String.Format("{0}", ts.Milliseconds);
                Console.WriteLine("RunTime " + elapsedTime);

                //вывод решения x
                for (int i = 0; i < X.N; i++) Console.WriteLine("X[{0}] = {1}", i + 1, X.Elem[i]);

                //подсчет и вывод погрешности
                var diff = X - X_true;
                var module1 = diff.modulus();
                var module2 = X_true.modulus();
                double inaccuracy = module1/module2;
                Console.WriteLine("Погрешность: "+inaccuracy);

            }
            catch (Exception E)
            {
                Console.WriteLine("\n*** Error! ***");
                Console.WriteLine("Method:  {0}",   E.TargetSite);
                Console.WriteLine("Message: {0}\n", E.Message);
            }
        }
    }
}