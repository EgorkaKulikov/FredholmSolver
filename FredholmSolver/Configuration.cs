using System;

namespace FredholmSolver
{
  class Configuration
  {
    //Количество узлов сетки
    public const int GridPoints = 15;

    //Отрезок аппроксимации (начало и конец сетки) и сдвиг для фиктивных точек
    public const double Left = 0.0;
    public const double Right = 1.0 * Math.PI / 2; 
    public static readonly double Eps = Math.Pow(10, -3);

    //Параметр регуляризации
    public static readonly double Alpha = Math.Pow(10, -10);

    //Функции, определяющие уравнение Фредгольма
    public static Func<double, double, double> K = (s, t) => Math.Sin(s) * Math.Cos(t);
    public static Func<double, double> F = Math.Sin;
    
    //Функция точного решения уравнения Фредгольма
    public static Func<double, double> U = s => 2 * Math.Sin(s);

    //Используемый метод аппроксимации
    public const ApproximationType ApproxType = ApproximationType.Sablonniere;
  }
}
