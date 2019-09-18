import pandas as pd

import unittest

from Phosphorpy.external import css


def test_download_light_curves():
    coord = pd.DataFrame(
        {
            'ra': [18.15629, 166.12397, 260.78809],
            'dec': [33.49162, 8.6418, 48.31078]
        }
    )
    lc = css.download_light_curves(coord['ra'].values,
                                   coord['dec'].values)
    # print(lc)
    lc = css.download_light_curves(coord['ra'],
                                   coord['dec'])
    # print(lc)


def test_download_light_curve():
    lc = css.download_light_curve(18.15629, 33.49162)
    print(lc)


def test_daily_average():

    lc = css.download_light_curve(18.15629, 33.49162)

    avg = css.daily_average(lc)

    assert len(avg) <= len(lc)


if __name__ == '__main__':
    unittest.main()
