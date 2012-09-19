// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "accuracy_options.hpp"

#include <stdexcept>
#include <algorithm>

namespace Fiber
{

namespace
{

struct LessOrEqual
{
    bool operator()(const std::pair<double, QuadratureOptions>& first,
                    const std::pair<double, QuadratureOptions>& second) const
    {
        return first.first <= second.first;
    }
};

struct Equal
{
    bool operator()(const std::pair<double, QuadratureOptions>& first,
                    const std::pair<double, QuadratureOptions>& second) const
    {
        return first.first == second.first;
    }
};

} // namespace

AccuracyOptionsEx::AccuracyOptionsEx()
{
    m_doubleRegular.push_back(std::make_pair(std::numeric_limits<double>::infinity(),
                                             QuadratureOptions()));
}

AccuracyOptionsEx::AccuracyOptionsEx(const AccuracyOptions& oldStyleOpts)
{
    m_singleRegular = oldStyleOpts.singleRegular;
    m_doubleRegular.push_back(std::make_pair(std::numeric_limits<double>::infinity(),
                                             oldStyleOpts.doubleRegular));
    m_doubleSingular = oldStyleOpts.doubleSingular;
}

const QuadratureOptions& AccuracyOptionsEx::singleRegular() const
{
    return m_singleRegular;
}

void AccuracyOptionsEx::setSingleRegular(
        int accuracyOrder, bool relativeToDefault)
{
    if (relativeToDefault)
        m_singleRegular.setRelativeQuadratureOrder(accuracyOrder);
    else
        m_singleRegular.setAbsoluteQuadratureOrder(accuracyOrder);
}

const QuadratureOptions& AccuracyOptionsEx::doubleRegular(
        double relativeDistance) const
{
    for (int i = m_doubleRegular.size() - 1; i >= 0; --i)
        if (relativeDistance <= m_doubleRegular[i].first)
            return m_doubleRegular[i].second;
    // should never happen
    throw std::runtime_error("AccuracyOptions::doubleRegular(): "
                             "internal error");
}

void AccuracyOptionsEx::setDoubleRegular(
        int accuracyOrder, bool relativeToDefault)
{
    m_doubleRegular.clear();
    m_doubleRegular.push_back(
                std::make_pair(std::numeric_limits<double>::infinity(),
                               QuadratureOptions(accuracyOrder, relativeToDefault)));
}

void AccuracyOptionsEx::setDoubleRegular(
        double maxNormalizedDistance1, int accuracyOrder1,
        int accuracyOrder2, bool relativeToDefault)
{
    m_doubleRegular.clear();
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance1,
                               QuadratureOptions(accuracyOrder1, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(std::numeric_limits<double>::infinity(),
                               QuadratureOptions(accuracyOrder2, relativeToDefault)));
    std::sort(m_doubleRegular.begin(), m_doubleRegular.end(), LessOrEqual());
    std::unique(m_doubleRegular.begin(), m_doubleRegular.end(), Equal());
}

void AccuracyOptionsEx::setDoubleRegular(
        double maxNormalizedDistance1, int accuracyOrder1,
        double maxNormalizedDistance2, int accuracyOrder2,
        int accuracyOrder3, bool relativeToDefault)
{
    m_doubleRegular.clear();
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance1,
                               QuadratureOptions(accuracyOrder1, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance2,
                               QuadratureOptions(accuracyOrder2, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(std::numeric_limits<double>::infinity(),
                               QuadratureOptions(accuracyOrder3, relativeToDefault)));
    std::sort(m_doubleRegular.begin(), m_doubleRegular.end(), LessOrEqual());
    std::unique(m_doubleRegular.begin(), m_doubleRegular.end(), Equal());
}

void AccuracyOptionsEx::setDoubleRegular(
        double maxNormalizedDistance1, int accuracyOrder1,
        double maxNormalizedDistance2, int accuracyOrder2,
        double maxNormalizedDistance3, int accuracyOrder3,
        int accuracyOrder4, bool relativeToDefault)
{
    m_doubleRegular.clear();
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance1,
                               QuadratureOptions(accuracyOrder1, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance2,
                               QuadratureOptions(accuracyOrder2, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance3,
                               QuadratureOptions(accuracyOrder3, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(std::numeric_limits<double>::infinity(),
                               QuadratureOptions(accuracyOrder4, relativeToDefault)));
    std::sort(m_doubleRegular.begin(), m_doubleRegular.end(), LessOrEqual());
    std::unique(m_doubleRegular.begin(), m_doubleRegular.end(), Equal());
}

void AccuracyOptionsEx::setDoubleRegular(
        double maxNormalizedDistance1, int accuracyOrder1,
        double maxNormalizedDistance2, int accuracyOrder2,
        double maxNormalizedDistance3, int accuracyOrder3,
        double maxNormalizedDistance4, int accuracyOrder4,
        int accuracyOrder5, bool relativeToDefault)
{
    m_doubleRegular.clear();
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance1,
                               QuadratureOptions(accuracyOrder1, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance2,
                               QuadratureOptions(accuracyOrder2, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance3,
                               QuadratureOptions(accuracyOrder3, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(maxNormalizedDistance4,
                               QuadratureOptions(accuracyOrder4, relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(std::numeric_limits<double>::infinity(),
                               QuadratureOptions(accuracyOrder5, relativeToDefault)));
    std::sort(m_doubleRegular.begin(), m_doubleRegular.end(), LessOrEqual());
    std::unique(m_doubleRegular.begin(), m_doubleRegular.end(), Equal());
}

void AccuracyOptionsEx::setDoubleRegular(
        const std::vector<double>& maxNormalizedDistances,
        const std::vector<int>& accuracyOrders,
        bool relativeToDefault)
{
    if (maxNormalizedDistances.size() != accuracyOrders.size() - 1)
        throw std::invalid_argument("AccuracyOptionsEx::setDoubleRegular(): "
                                    "maxNormalizedDistances must have one "
                                    "element less than accuracyOrders");
    m_doubleRegular.clear();
    for (size_t i = 0; i < maxNormalizedDistances.size(); ++i)
        m_doubleRegular.push_back(
                    std::make_pair(maxNormalizedDistances[i],
                                   QuadratureOptions(accuracyOrders[i],
                                                     relativeToDefault)));
    m_doubleRegular.push_back(
                std::make_pair(std::numeric_limits<double>::infinity(),
                               QuadratureOptions(accuracyOrders.back(),
                                                 relativeToDefault)));
    std::sort(m_doubleRegular.begin(), m_doubleRegular.end(), LessOrEqual());
    std::unique(m_doubleRegular.begin(), m_doubleRegular.end(), Equal());
}

const QuadratureOptions& AccuracyOptionsEx::doubleSingular() const
{
    return m_doubleSingular;
}

void AccuracyOptionsEx::setDoubleSingular(
        int accuracyOrder, bool relativeToDefault)
{
    if (relativeToDefault)
        m_doubleSingular.setRelativeQuadratureOrder(accuracyOrder);
    else
        m_doubleSingular.setAbsoluteQuadratureOrder(accuracyOrder);
}

} // namespace Fiber
