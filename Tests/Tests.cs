namespace Magnetostatics.Tests;

public interface ITest
{
    double Az(Point2D point);

    double J(Point2D point);
}

public class Test1 : ITest
{
    public double Az(Point2D point) => point.X + point.Y;

    public double J(Point2D point) => 0.0;
}