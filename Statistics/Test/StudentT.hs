{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DeriveDataTypeable #-}
-- |
-- Module    : Statistics.Test.StudentT
-- Copyright : (c) 2013 Chris Taylor
-- License   : BSD3
--
-- Maintainer  : bos@serpentine.com
-- Stability   : experimental
-- Portability : portable
--
-- There are several variations on Student't t-test. They are all parametric
-- tests. The one-sample t-test is used to determine whether a sample has a
-- nonzero mean. The two-sample t-test is used to determine whether two independent
-- samples have different means. The paired t-test is used to determine whether
-- two dependent samples have different means.
module Statistics.Test.StudentT (
    ttestOneSample
  , ttestTwoSample
  , ttestPaired
    -- * Data types
  , VarianceAssumption(..)
  , TestType(..)
  , TestResult(..)
  ) where

import qualified Data.Vector.Generic as G
import           Data.Typeable (Typeable)

import           Statistics.Distribution
import           Statistics.Distribution.StudentT
import           Statistics.Test.Types
import qualified Statistics.Sample as Sample

-- |There are two versions of the two-sample t-test, depending on whether
-- the samples can be assumed to have equal variance or not. This allows
-- you to switch between the two versions.
data VarianceAssumption = EqualVariance
                        | UnequalVariance
                          deriving (Eq,Ord,Show,Typeable)

-- | Computes the sample p-value for the Student's t distribution with n
-- degrees of freedom. The computed value can be either one-tailed or
-- two-tailed.
computePval :: TestType -> Double -> Double -> Double
computePval testType dof t = case testType of
  OneTailed -> complCumulative d t
  TwoTailed -> complCumulative d t' + cumulative d (negate t')
 where
  d  = studentT dof
  t' = abs t

-- | Performs a t-test of the hypothesis that the data in a vector come
-- from a distribution with mean zero.
ttestOneSample :: (G.Vector v Double)
  => TestType     -- ^ One or two-tailed test?
  -> Double       -- ^ Significance level (e.g. 0.05)
  -> v Double     -- ^ Vector of observations
  -> TestResult
ttestOneSample testType p vec
  | n < 1          = error $ "Statistics.Test.StudentT.ttestOneSample: too short data sample"
  | p > 0 && p < 1 = significant $ computePval testType n t < p
  | otherwise      = error $ "Statistics.Test.StudentT.ttestOneSample: bad p-value: " ++ show p
  where
    n = fromIntegral $ G.length vec
    m = Sample.mean vec
    s = Sample.stdDev vec
    t = m / s * sqrt n

-- | Performs a paired t-test of the hypothesis that two matched samples
-- come from distributions with equal means. The difference between each
-- pair of observation is assumed to come from a normal distribution with
-- unknown variance. The two samples must have the same length.
ttestPaired :: (G.Vector v Double)
  => TestType     -- ^ One or two-tailed test?
  -> Double       -- ^ Significance level (e.g. 0.05)
  -> v Double     -- ^ First vector of observations
  -> v Double     -- ^ Second vector of observations
  -> TestResult
ttestPaired testType p x y
  | G.length x /= G.length y = error $ "Statistics.Test.StudentT.ttestPaired: unbalanced data samples"
  | n < 1                    = error $ "Statistics.Test.StudentT.ttestPaired: too short data sample"
  | p > 0 && p < 1           = significant $ computePval testType n t < p
  | otherwise                = error $ "Statistics.Test.StudentT.ttestPaired: bad p-value: " ++ show p
  where
    v = G.zipWith (-) x y
    n = fromIntegral $ G.length v
    m = Sample.mean v
    s = Sample.stdDev v
    t = m / s * sqrt n

ttestTwoSample :: (G.Vector v Double)
  => TestType             -- ^ One or two-tailed test?
  -> VarianceAssumption   -- ^ Assume equal variances?
  -> Double               -- ^ Significance level (e.g. 0.05)
  -> v Double             -- ^ First vector of observations
  -> v Double             -- ^ Second vector of observations
  -> TestResult
ttestTwoSample testType varAssumption p x y
  | n < 1 || m < 1 = error $ "Statistics.Test.StudentT.ttestTwoSample: too short data sample"
  | p > 0 && p < 1 = significant $ computePval testType dof t < p
  | otherwise      = error $ "Statistics.Test.StudentT.ttestTwoSample: bad p-value: " ++ show p
  where
    n = fromIntegral $ G.length x
    m = fromIntegral $ G.length y
    mx = Sample.mean x
    my = Sample.mean y
    sx = Sample.varianceUnbiased x
    sy = Sample.varianceUnbiased y
    dof = case varAssumption of
      EqualVariance   -> n + m - 2
      UnequalVariance -> (sx/n + sy/m)^2 / ((sx/n)^2/(n-1) + (sy/m)^2/(m-1))
    denom = case varAssumption of
      EqualVariance   -> sqrt $ ((n-1) * sx + (m-1) * sy) / (n + m - 2) * (1/n + 1/m)
      UnequalVariance -> sqrt $ (sx/n + sy/n)
    t = (mx - my) / denom



