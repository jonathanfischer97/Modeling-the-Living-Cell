function calc_gauss(mu,sigma)
    xmin = mu-(sigma*4)
    xmax = mu+(sigma*4)
    nPts = 100
    dx = (xmax-xmin)/100
    xValues = collect(xmin:dx:xmax)

    f(x) = 1/sqrt(2*pi*sigma^2)*exp(-1/(2*sigma^2)*(x-mu)^2)
    yValues = map(f, xValues)

    (xValues, yValues)
end

