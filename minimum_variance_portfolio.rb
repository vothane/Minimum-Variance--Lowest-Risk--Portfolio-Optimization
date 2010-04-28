require 'matrix'

class MinimumVariancePortfolio

  attr_accessor :expected_mean_of_returns, :variance_of_returns, :weights

  def initialize( expected_mean_of_returns, variance_of_returns )
    if (expected_mean_of_returns.instance_of? Matrix) and (variance_of_returns.instance_of? Matrix)
      @expected_mean_of_returns = expected_mean_of_returns
      @variance_of_returns      = variance_of_returns
    elsif (expected_mean_of_returns.instance_of? Array) and (variance_of_returns.instance_of? Array)
      @expected_mean_of_returns = Matrix[expected_mean_of_returns]
      @variance_of_returns      = Matrix.diagonal( *variance_of_returns )
    else
      raise ArgumentError, "Wrong Arguments"
    end
  end

  def calculate( return_level )
    raise ArgumentError, "The target return level must be between 0 and 1" if return_level > 1.0
    
    mean     = @expected_mean_of_returns.transpose
    variance = @variance_of_returns
    n        = mean.to_a.length

    a = (ones( 1, n ) * variance.inverse * mean).to_a.flatten.first
    b = (mean.transpose * variance.inverse * mean).to_a.flatten.first
    c = (ones( 1, n ) * variance.inverse * ones( n, 1 )).to_a.flatten.first
    d = (Matrix[[b,a],[a,c]]).determinant
    e = (1/d) * (b * ones( 1, n ) - a * mean.transpose) * variance.inverse
    f = (1/d) * (c * mean.transpose - a * ones( 1, n )) * variance.inverse
  
    allocated_weights = (e + (f * return_level))
    @weights = allocated_weights.to_a.flatten
    return minimum_variance_portfolio_exists?
  end

  private

  def minimum_variance_portfolio_exists?
    @weights.min > 0.0
  end

  def ones( rows, columns )
    a = []
    rows.times { a << Array.new( columns, 1 ) }
    Matrix[*a]
  end
end