using System;

namespace Com_Methods
{
    /// <summary>
    /// преобразования Хаусхолдера
    /// </summary>
    class Householder_Transformation
    {
        /// <summary>
        /// реализация процедуры ортогонализации по методу отражений Хаусхолдера
        /// </summary>
        /// <param name="A - исходная матрица"></param>
        /// <param name="Q - ортогональная матрица преобразований"></param>
        /// <param name="R - результат"></param>
        public static void Householder_Orthogonalization(Matrix A, Matrix Q, Matrix R)
        {
            //инициализация вектора отражения
            Vector p = new Vector(A.M);

            //вспомогательные переменные
            double s, beta, mu;

            //запись матрицы А в R
            R.Copy(A);

            //алгоритм отражений Хаусхолдера
            for (int i = 0; i < R.N - 1; i++)
            {
                //находим квадрат нормы столбца для обнуления
                s = 0;
                for (int I = i; I < R.M; I++) s += Math.Pow(R.Elem[I][i], 2);

                //если есть ненулевые элементы под диагональю, тогда 
                //квадрат нормы вектора для обнуления не совпадает с квадратом диагонального элемента
                if (Math.Sqrt(Math.Abs(s - R.Elem[i][i] * R.Elem[i][i])) > CONST.EPS)
                {
                    //выбор знака слагаемого beta = sign(-x1)
                    if (R.Elem[i][i] < 0) beta = Math.Sqrt(s);
                    else beta = -Math.Sqrt(s);

                    //вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
                    mu = 1.0 / beta / (beta - R.Elem[i][i]);

                    //формируем вектор p
                    for (int I = 0; I < R.M; I++) { p.Elem[I] = 0; if (I >= i) p.Elem[I] = R.Elem[I][i]; }

                    //изменяем диагональный элемент
                    p.Elem[i] -= beta;

                    //вычисляем новые компоненты матрицы A = Hm * H(m-1) ... * A
                    for (int m = i; m < R.N; m++)
                    {
                        //произведение S = At * p
                        s = 0;
                        for (int n = i; n < R.M; n++) { s += R.Elem[n][m] * p.Elem[n]; }
                        s *= mu;
                        //A = A - 2 * p * (At * p)^t / ||p||^2
                        for (int n = i; n < R.M; n++) { R.Elem[n][m] -= s * p.Elem[n]; }
                    }

                    //вычисляем новые компоненты матрицы Q = Q * H1 * H2 * ...
                    for (int m = 0; m < R.M; m++)
                    {
                        //произведение Q * p
                        s = 0;
                        for (int n = i; n < R.M; n++) { s += Q.Elem[m][n] * p.Elem[n]; }
                        s *= mu;
                        //Q = Q - p * (Q * p)^t
                        for (int n = i; n < R.M; n++) { Q.Elem[m][n] -= s * p.Elem[n]; }
                    }
                }
            }
        }
    }
}
