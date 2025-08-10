using System;

namespace FredholmSolver
{
  class Program
  {
    public static void Main(string[] args)
    {
      Program program = new Program();

      double error = program.CreateExperiment();

      Console.WriteLine("-------------------------------------------------");
      Console.WriteLine(error);
      Console.ReadKey();
    }

    /// <summary>
    /// Вычисление ошибки аппроксимации: строим значения функции на сетке, 
    /// в десять раз меньшей исходной, и сравниваем с реальными значениями функции.
    /// </summary>
    private double CreateExperiment()
    {
      ApproximationEngine engine = new ApproximationEngine();
      
      //шаг мелкой сетки
      double step = 1.0 * (Configuration.Right - Configuration.Left) / Configuration.GridPoints / 10;

      double error = 0.0;
      Console.WriteLine("Point  Value   Appr    ApprSloan   Error    SloanError");
      Console.WriteLine("---------------------------------------------------------");

      for (double point = Configuration.Left; point < Configuration.Right; point += step)
      {
        //вычисляем реальное значение, приближённое, а также приближённое, уточнённое по Слоану
        double realValue = Configuration.U(point);
        double approxValue = engine.ApproxInPoint(point);
        double approxSloanValue = engine.ApproxSloanInPoint(point);

        //вычисляем ошибки методов приближения в данной точке
        double pointError = Math.Abs(realValue - approxValue);
        double pointSloanError = Math.Abs(realValue - approxSloanValue);

        Console.WriteLine(
          "{0:0.00}  {1:0.00000} {2:0.00000}  {3:0.00000}    {4:0.00000}    {5:0.00000}",
          point, realValue, approxValue, approxSloanValue, pointError, pointSloanError);

        //вычисляем максимальное значение ошибки
        if (pointSloanError > error)
        {
          error = pointSloanError;
        }
      }

      return error;
    }
  }
}
