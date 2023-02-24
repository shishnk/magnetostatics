namespace Magnetostatics.Tests;

public interface ITest
{
    double Az(Point2D point);

    double J(Point2D point);
}