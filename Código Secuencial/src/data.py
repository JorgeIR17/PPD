from typing import Iterable, Tuple, Union

import numpy as np
import pandas as pd

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


def crop_dataframe(df: pd.DataFrame, bounding_box: Iterable[float]) -> pd.DataFrame:
    """
    Crop a DataFrame to a bounding box.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to be cropped.
    bounding_box : Iterable[float]
        Bounding box of the DataFrame, defined by (latmin, lonmin, latmax, lonmax).

    Returns
    -------
    pd.DataFrame
        Cropped DataFrame.
    """
    latmin, lonmin, latmax, lonmax = bounding_box
    return df.loc[
        (df.index > latmin) & (df.index < latmax),
        (df.columns > lonmin) & (df.columns < lonmax),
    ]


class Ocean:
    def __init__(
        self, df_vo: Union[str, pd.DataFrame], df_uo: Union[str, pd.DataFrame]
    ):
        """
        Initialize the Ocean object.

        Parameters
        ----------
        df_vo : Union[str, pd.DataFrame]
            DataFrame containing the v (latitude) component of the ocean currents.
            If a string is passed, it will be consider as the path to a CSV file.
        df_uo : Union[str, pd.DataFrame]
            DataFrame containing the u (longitude) component of the ocean currents.
            If a string is passed, it will be consider as the path to a CSV file.
        """

        if isinstance(df_uo, str):
            df_uo = pd.read_csv(df_uo, index_col=0)
            df_uo.columns = df_uo.columns.astype(float)
        if isinstance(df_vo, str):
            df_vo = pd.read_csv(df_vo, index_col=0)
            df_vo.columns = df_vo.columns.astype(float)

        self.uo = df_uo
        self.vo = df_vo
        self.mod = pd.DataFrame(
            np.sqrt(df_uo.to_numpy() ** 2 + df_vo.to_numpy() ** 2),
            columns=df_uo.columns,
            index=df_uo.index,
        )

    def crop_data(self, bounding_box: Iterable[float]):
        """
        Crop the ocean data to a bounding box.

        Parameters
        ----------
        bounding_box : Iterable[float]
            Bounding box of the DataFrame, defined by (latmin, lonmin, latmax, lonmax).
        """

        self.uo = crop_dataframe(self.uo, bounding_box)
        self.vo = crop_dataframe(self.vo, bounding_box)
        self.mod = np.sqrt(self.uo**2 + self.vo**2)

    def get_currents(self, lat: np.array, lon: np.array) -> Tuple[np.array]:
        """
        Get the currents at the given latatitude, longitude coordinate.

        Parameters
        ----------
        lat : np.array
            Array containing the latitude of the points.
        lon : np.array
            Array containing the longitude of the points.

        Returns
        -------
        Tuple[np.array]
            Tuple containing the v and u components of the
        """

        y = np.abs(self.uo.index.values - lat).argmin()
        x = np.abs(self.uo.columns.values - lon).argmin()

        uo = self.uo.iloc[y, x]
        vo = self.vo.iloc[y, x]

        return vo, uo

    def plot_currents(
        self,
        cmap: str = "viridis",
        land_color: str = "0.8",
        water_color_background: bool = False,
        data_on_top: bool = False,
        add_cb: bool = True,
        draw_coastlines: bool = False,
        vmin: float = None,
        vmax: float = None,
        cb_title: str = "Currents speed (m/s)",
        fig=None,
        ax=None,
    ) -> tuple:
        """
        Plot a map using matplotlib.

        Parameters
        ----------
        cmap : str, optional
            Colormap to use, by default "viridis".
        add_cb : bool, optional
            If True, a colorbar will be added, by default True.
        draw_coastlines : bool, optional
            If True, coastlines will be drawn, by default False.
        vmin: float, optional
            Minimum value of the colorbar, by default None.
        vmax: float, optional
            Maximum value of the colorbar, by default None.
        land_color : str, optional
            Color of the inland, by default "0.8".
            "0.8" is a light gray.
        water_color_background: bool, optional
            If True the water will be colored as the minimym value color of the cmap, by default False.
            Else, the water will be colored as the land.
        data_on_top : bool, optional
            If True, the data will be plotted on top of the map, by default False.
        cb_title : str, optional
            Title of the colorbar, by default "Currents (m/s)".
            If it is None, the colorbar will not have a title.
        fig : plt.figure, optional
            Figure where the map will be drawn, by default None.
        ax : Axes, optional
            Axes where the map will be drawn, by default None.

        Returns
        -------
        tuple
            Basemap object, figure and axis of the map.
        """
        latmin, lonmin = self.uo.index.min(), self.uo.columns.min()
        latmax, lonmax = self.uo.index.max(), self.uo.columns.max()

        m = Basemap(
            llcrnrlat=latmin,
            urcrnrlat=latmax,
            llcrnrlon=lonmin,
            urcrnrlon=lonmax,
            resolution="h",
            ax=ax,
        )

        # TODO: use pcolormesh instead of imshow.
        # This should allow to make the bounding_box bigger than the data
        data_mod = np.flip(self.mod.to_numpy(), axis=0)
        data_mod = np.nan_to_num(data_mod, 0)
        cs = m.imshow(
            data_mod,
            cmap=cmap,
            zorder=20 if data_on_top else 0,
            vmin=vmin,
            vmax=vmax,
        )

        if add_cb:
            cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
            water_color = (
                plt.cm.get_cmap(cmap, 20)(0) if water_color_background else land_color
            )
            if cb_title is not None:
                cb.set_label(cb_title)
        else:
            water_color = land_color

        m.drawlsmask(
            land_color=land_color, ocean_color=(0, 0, 0, 0), lakes=True, resolution="h"
        )
        m.fillcontinents(color=land_color, lake_color=water_color, zorder=10)

        if draw_coastlines:
            m.drawcoastlines(zorder=30)

        return m, fig, ax
