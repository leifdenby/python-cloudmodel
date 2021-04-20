import matplotlib.pyplot as plot
import matplotlib.gridspec as gridspec
import numpy as np
import warnings
import math

try:
    from tephigram_python.tephigram_plotter import Tephigram

except ImportError:
    Tephigram = None

from ..reference import parameterisations
from .. import Var


def profile_plot(F, z, Te):
    r = F[:, Var.r]
    w = F[:, Var.w]
    T = F[:, Var.T]

    plot.figure(figsize=(16, 6))
    plot.subplot(131)
    plot.plot(r, z, marker="x")
    plot.ylabel("height [m]")
    plot.xlabel("radius [m]")
    plot.grid(True)
    plot.xlim(0.0, None)
    plot.subplot(132)
    plot.plot(T, z, marker="x", label="in-cloud")
    plot.plot(Te(z), z, marker="x", label="environment")
    plot.legend()
    plot.xlabel("temperature [K]")
    plot.ylabel("height [m]")
    plot.grid(True)
    plot.subplot(133)
    plot.plot(w, z, marker="x")
    plot.ylabel("height [m]")
    plot.xlabel("vertical velocity [m/s]")
    plot.grid(True)

    return plot


def hydrometeor_profile_plot(F, z, Te, p_e):
    p = p_e
    T = F[:, Var.T]
    q_v = F[:, Var.q_v]
    q_l = F[:, Var.q_l]
    q_r = F[:, Var.q_r]

    T = [float(t_) for t_ in T]

    q_v__sat = parameterisations.qv_sat(T=T, p=p)

    plot.figure(figsize=(16, 6))
    plot.subplot(131)
    plot.plot(q_v, z, marker="x", label="in-cloud")
    plot.plot(q_v__sat, z, marker="x", label="saturation (in-cloud temp)")
    plot.ylabel("height [m]")
    plot.xlabel("water vapor specific mass [kg/kg]")
    plot.legend()
    plot.grid(True)
    plot.subplot(132)
    plot.plot(q_l, z, marker="x", label="in-cloud")
    plot.ylabel("height [m]")
    plot.xlabel("liquid water specific mass [kg/kg]")
    plot.grid(True)
    plot.subplot(133)
    plot.plot(q_r, z, marker="x", label="in-cloud")
    plot.ylabel("height [m]")
    plot.xlabel("rain water specific mass [kg/kg]")
    plot.grid(True)


def _get_var(v, profile):
    z = profile.z

    ref_plot_func = None
    if v in Var.names:
        profile_data = profile.F[:, Var.names.index(v)]
    elif v in ["d_lse", "lse"]:
        constants = profile.cloud_model.constants
        cp_d = constants.cp_d
        cp_v = constants.cp_v
        L_v = constants.L_v
        L_s = constants.L_s
        g = constants.g

        def lse_env_f(z_):
            p = profile.cloud_model.environment.p(z_)
            T_e = profile.cloud_model.environment.temp(z_)

            qv_sat__f = parameterisations.ParametersationsWithSpecificConstants(
                constants=constants
            ).pv_sat.qv_sat
            qv_e__sat = qv_sat__f(T=T_e, p=p)
            rh_e = profile.cloud_model.environment.rel_humidity(z_)
            qv_e = rh_e * qv_e__sat
            qd_e = 1.0 - qv_e
            ql_e, qi_e = 0.0, 0.0

            c_em_p = cp_d * qd_e + cp_v * (qv_e + ql_e + qi_e)

            lse_e = c_em_p * T_e - ql_e * L_v + g * z_

            return lse_e

        T_c = profile.F[:, Var.T]
        qv_c = profile.F[:, Var.q_v]
        qr_c = profile.F[:, Var.q_r]
        ql_c = profile.F[:, Var.q_l]
        qi_c = profile.F[:, Var.q_i]
        qd_c = 1.0 - qv_c - ql_c - qi_c

        c_cm_p = cp_d * qd_c + cp_v * (qv_c + ql_c + qi_c)
        lse_c = c_cm_p * T_c - ql_c * L_v + z * g

        if v == "d_lse":
            z = profile.z
            lse_e = lse_env_f(z)

            Ds = lse_c - lse_e
            profile_data = Ds
        elif v == "lse":
            z_ = np.linspace(0.0, profile.z.max(), 100)
            lse_e = lse_env_f(z_)
            profile_data = lse_c
            ref_plot_func = lambda **kwargs: plot.plot(
                lse_e, z_, marker="", label="environment", **kwargs
            )
    elif v in ["d_mse", "mse"]:
        constants = profile.cloud_model.constants
        cp_d = constants.cp_d
        cp_l = constants.cp_l
        L_v = constants.L_v
        g = constants.g

        def mse_env_f(z_):
            p = profile.cloud_model.environment.p(z_)
            T_e = profile.cloud_model.environment.temp(z_)

            qv_sat__f = parameterisations.ParametersationsWithSpecificConstants(
                constants=constants
            ).pv_sat.qv_sat
            qv_e__sat = qv_sat__f(T=T_e, p=p)
            rh_e = profile.cloud_model.environment.rel_humidity(z_)
            qv_e = rh_e * qv_e__sat
            qd_e = 1.0 - qv_e
            ql_e, qi_e = 0.0, 0.0

            c_em_p = cp_d * qd_e + cp_l * (qv_e + ql_e + qi_e)

            mse_e = c_em_p * T_e + qv_e * L_v + g * z_

            return mse_e

        T_c = profile.F[:, Var.T]
        qv_c = profile.F[:, Var.q_v]
        qr_c = profile.F[:, Var.q_r]
        ql_c = profile.F[:, Var.q_l]
        qi_c = profile.F[:, Var.q_i]
        qd_c = 1.0 - qv_c - ql_c - qi_c - qr_c

        c_cm_p = cp_d * qd_c + cp_l * (qv_c + ql_c + qr_c + qi_c)
        mse_c = c_cm_p * T_c + qv_c * L_v + z * g

        if v == "d_mse":
            z = profile.z
            mse_e = mse_env_f(z)

            Ds = mse_c - mse_e
            profile_data = Ds
        elif v == "mse":
            z_ = np.linspace(0.0, profile.z.max(), 100)
            mse_e = mse_env_f(z_)
            profile_data = mse_c
            ref_plot_func = lambda **kwargs: plot.plot(
                mse_e / 1000.0, z_, marker="", label="environment", **kwargs
            )
    elif v == "Sw":
        profile_data = None
    elif v == "rho_c":
        z = profile.z
        p = profile.cloud_model.environment.p(z)
        T = profile.F[:, Var.T]
        qv_c = profile.F[:, Var.q_v]
        qr_c = profile.F[:, Var.q_r]
        ql_c = profile.F[:, Var.q_l]
        qi_c = profile.F[:, Var.q_i]
        qd_c = 1.0 - qv_c - ql_c - qi_c
        profile_data = profile.cloud_model.cloud_mixture_density(
            p=p, T_c=T, qd_c=qd_c, qv_c=qv_c, ql_c=ql_c, qi_c=qi_c, qr_c=qr_c
        )
        ref_plot_func = lambda **kwargs: plot.plot(
            profile.cloud_model.environment.rho(profile.z),
            profile.z,
            marker="",
            label="environment",
            **kwargs
        )
    elif v == "d_rho":
        z = profile.z
        p = profile.cloud_model.environment.p(z)
        T_c = profile.F[:, Var.T]
        qv_c = profile.F[:, Var.q_v]
        qr_c = profile.F[:, Var.q_r]
        ql_c = profile.F[:, Var.q_l]
        qi_c = profile.F[:, Var.q_i]
        qd_c = 1.0 - qv_c - ql_c - qi_c - qr_c
        rho_c = profile.cloud_model.cloud_mixture_density(
            p=p, T_c=T_c, qd_c=qd_c, qv_c=qv_c, ql_c=ql_c, qi_c=qi_c, qr_c=qr_c
        )
        rho_e = profile.cloud_model.environment.rho(profile.z)

        profile_data = rho_c - rho_e
    elif v == "d_qv":
        T = profile.F[:, Var.T]
        p = profile.F[:, Var.p]
        qv_c = profile.F[:, Var.q_v]
        z = profile.z
        p_e = profile.cloud_model.environment.p(z)
        constants = profile.cloud_model.constants

        T_e = profile.cloud_model.environment.temp(z)
        qv_e__sat = parameterisations.ParametersationsWithSpecificConstants(
            constants=constants
        ).pv_sat.qv_sat(T=T_e, p=p)
        rh_e = profile.cloud_model.environment.rel_humidity(z)
        qv_e = rh_e * qv_e__sat

        profile_data = (qv_c - qv_e) * 1000

    elif v == "RH":
        # handled further below
        profile_data = None

    elif hasattr(profile, "extra_vars") and v in profile.extra_vars:
        profile_data = profile.extra_vars[v]
        if len(profile_data) == len(profile.z):
            pass
        else:
            try:
                z = profile.extra_vars["t_substeps"]
            except KeyError:
                warnings.warn(
                    "Had to skip plotting `{}` because too many datapoints were"
                    " found, probably using a sub-stepping integration method".format(v)
                )

        if v == "r_c":
            profile_data = 1.0e6 * np.array(profile_data)
    else:
        raise NotImplementedError(v)

    return profile_data, ref_plot_func


def plot_profiles(
    profiles,
    variables=["r", "w", "T", "q_v", "q_l"],
    initial_condition=None,
    labels_ncol=1,
    label_f=None,
    col_max=3,
    fig=None,
):
    n = len(variables)
    c = n > col_max and col_max or n
    r = int(math.ceil(float(n) / c))

    gs = gridspec.GridSpec(r, c)

    if fig is None:
        fig = plot.figure(figsize=(6 * c, 4 * r))

    lines = []
    for n, (v, s) in enumerate(zip(variables, list(gs))):

        tephigram = None
        if v == "T__tephigram":
            if Tephigram is None:
                raise Exception(
                    "Can't plot tephigram without `tephigram_python` installed"
                )
            tephigram = Tephigram(fig=fig, subplotshape=(r, c, n + 1))
            v = "T"
        else:
            plot.subplot(s)

        ref_lines = []

        scale_by_max = False
        d_max = 0.0

        for n_profile, profile in enumerate(profiles):
            z = profile.z

            extra_var = False

            if tephigram is not None:
                z = profile.z
                p = profile.cloud_model.environment.p(z)
                T = profile.F[:, Var.T]
                kwargs = {"P": p / 100.0, "T": T - 273.15, "marker": "."}
                if len(lines) > 0:
                    kwargs["color"] = lines[n_profile].get_color()
                tephigram.plot_temp(**kwargs)

                RH = profile.F[:, Var.q_v] / parameterisations.pv_sat.qv_sat(T=T, p=p)
                kwargs = {"P": p / 100.0, "T": T - 273.15, "RH": RH}
                (RH_line,) = tephigram.plot_RH(**kwargs)
                RH_line.set_marker("")
                RH_line.set_linestyle("-")

                def ref_plot_func():
                    lines = []
                    z_ = np.linspace(0.0, z.max(), 100)
                    p_ = profile.cloud_model.environment.p(z_)

                    T_e = profile.cloud_model.environment.temp(z_)
                    kwargs = {"P": p_ / 100.0, "T": T_e - 273.15, "marker": "."}
                    (T_line,) = tephigram.plot_temp(**kwargs)
                    T_line.set_linestyle("-")
                    T_line.set_marker("")
                    T_line.set_label("env. T")
                    lines.append(T_line)

                    RH = profile.cloud_model.environment.rel_humidity(z_)
                    kwargs = {"P": p_ / 100.0, "T": T_e - 273.15, "RH": RH}
                    (RH_line,) = tephigram.plot_RH(**kwargs)
                    RH_line.set_linestyle("-")
                    RH_line.set_marker("")
                    RH_line.set_label("env. RH")
                    lines.append(RH_line)

                    # kwargs = { 'P': p/100., 'T': T-273.15, 'marker': '', 'color': 'black', 'label': 'environment',
                    # 'with_height_markers': z, 'marker_interval': 500, }
                    # return tephigram.plot_temp(**kwargs)
                    return lines

                plot.title("Tephigram")

            else:
                if label_f is not None:
                    label = label_f(profile)
                else:
                    label = str(profile)

                profile_data, ref_plot_func = _get_var(v=v, profile=profile)

                if not profile_data is None:
                    if v in "q_v q_l q_r q_i q_pr d_rho".split():
                        profile_data = 1000.0 * np.array(profile_data)
                    elif v in "lse d_lse mse d_mse".split():
                        # J/kg => kJ/kg
                        profile_data = np.array(profile_data) / 1000.0

                    profile_line = plot.plot(
                        profile_data,
                        z,
                        label=label,
                        marker=".",
                        linestyle="",
                    )

                    if n == 0:
                        lines += profile_line

                    d_max = np.nanmax([np.nanmax(profile_data), d_max])

                plot.grid(True)

                if v == "T":
                    plot.xlabel("temperature [K]")
                    plot.ylabel("height [m]")
                    ref_plot_func = lambda **kwargs: plot.plot(
                        profile.cloud_model.environment.temp(profile.z),
                        profile.z,
                        marker="",
                        label="environment",
                        **kwargs
                    )
                elif v == "r":
                    plot.ylabel("height [m]")
                    plot.xlabel("radius [m]")
                    plot.xlim(0.0, None)
                    scale_by_max = True
                elif v == "w":
                    plot.ylabel("height [m]")
                    plot.xlabel("vertical velocity [m/s]")
                elif v == "q_v":
                    plot.ylabel("height [m]")
                    plot.xlabel("water vapor specific mass [g/kg]")
                    scale_by_max = True

                    T = profile.F[:, Var.T]
                    p = profile.F[:, Var.p]
                    z = profile.z
                    p_e = profile.cloud_model.environment.p(z)
                    constants = profile.cloud_model.constants
                    q_v__sat = parameterisations.ParametersationsWithSpecificConstants(
                        constants=constants
                    ).pv_sat.qv_sat(T=T, p=p)
                    color = lines[n_profile].get_color()
                    plot.plot(q_v__sat * 1000.0, z, marker="", color=color, label="")

                    T_e = profile.cloud_model.environment.temp(z)
                    qv_e__sat = parameterisations.ParametersationsWithSpecificConstants(
                        constants=constants
                    ).pv_sat.qv_sat(T=T_e, p=p)
                    if not hasattr(profile.cloud_model.environment, "rel_humidity"):
                        rh_e = np.zeros_like(T_e)
                    else:
                        rh_e = profile.cloud_model.environment.rel_humidity(z)
                    qv_e = rh_e * qv_e__sat
                    ref_plot_func = lambda **kwargs: plot.plot(
                        qv_e * 1000.0, z, marker="", label="environment", **kwargs
                    )
                elif v == "d_qv":
                    plot.ylabel("height [m]")
                    plot.xlabel(r"$\Delta$ water vapor specific mass [mg/kg]")
                    scale_by_max = True
                elif v == "Sw":
                    plot.ylabel("height [m]")
                    plot.xlabel("super saturation [%]")

                    T = profile.F[:, Var.T]
                    p = profile.F[:, Var.p]
                    z = profile.z
                    constants = profile.cloud_model.constants
                    q_v__sat = parameterisations.ParametersationsWithSpecificConstants(
                        constants=constants
                    ).pv_sat.qv_sat(T=T, p=p)
                    q_v = profile.F[:, Var.q_v]
                    color = lines[n_profile].get_color()
                    Sw = (q_v / q_v__sat - 1.0) * 100.0
                    plot.plot(Sw, z, marker=".", color=color, label="", linestyle="")
                    plot.xlim(-5, 5)
                elif v == "RH":
                    plot.ylabel("height [m]")
                    plot.xlabel("relative humidity [%]")

                    T = profile.F[:, Var.T]
                    p = profile.F[:, Var.p]
                    z = profile.z
                    constants = profile.cloud_model.constants
                    q_v__sat = parameterisations.ParametersationsWithSpecificConstants(
                        constants=constants
                    ).pv_sat.qv_sat(T=T, p=p)
                    q_v = profile.F[:, Var.q_v]
                    color = lines[n_profile].get_color()
                    RH = q_v / q_v__sat * 100.0
                    plot.plot(RH, z, marker=".", color=color, label="", linestyle="")

                    z_ = np.linspace(0.0, z.max(), 100)
                    if hasattr(profile.cloud_model.environment, "rel_humidity"):
                        RH_e = 100 * profile.cloud_model.environment.rel_humidity(z_)
                        plot.plot(
                            RH_e, z_, label="", marker="", color=color, linestyle="-"
                        )

                elif v == "q_l":
                    plot.ylabel("height [m]")
                    plot.xlabel("liquid water specific mass [g/kg]")
                    scale_by_max = True
                elif v == "q_r":
                    plot.ylabel("height [m]")
                    plot.xlabel("rain water specific mass [g/kg]")
                    scale_by_max = True
                elif v == "q_pr":
                    plot.ylabel("height [m]")
                    plot.xlabel("precipitated rain water [g/kg]")
                    scale_by_max = True
                elif v == "d_lse":
                    plot.ylabel("height [m]")
                    plot.xlabel(
                        "Liquid static energy difference to environment [kJ/kg]"
                    )
                elif v == "lse":
                    plot.ylabel("height [m]")
                    plot.xlabel("Liquid static energy [kJ/kg]")
                elif v == "d_mse":
                    plot.ylabel("height [m]")
                    plot.xlabel("Moist static energy difference to environment [kJ/kg]")
                elif v == "mse":
                    plot.ylabel("height [m]")
                    plot.xlabel("Moist static energy [kJ/kg]")
                elif v == "r_c":
                    plot.ylabel(r"height [$m$]")
                    plot.xlabel("cloud-droplet radius [$\mu m$]")
                elif v == "Nc":
                    plot.ylabel("height [m]")
                    plot.xlabel("cloud droplet number [1/m^3]")
                elif v == "rho_c":
                    plot.ylabel("height [m]")
                    plot.xlabel("in-cloud density [kg/m3]")
                elif v == "d_rho":
                    plot.ylabel("height [m]")
                    plot.xlabel("density difference to environment [g/m^3]")
                elif extra_var:
                    plot.ylabel("height [m]")
                    plot.xlabel("{} [?]".format(v))
                else:
                    raise NotImplementedError

        if scale_by_max and d_max != 0.0:
            dx = 0.1 * d_max
            plot.xlim(0.0 - dx, d_max + dx)
        ylim = plot.gca().get_ylim()
        if ylim[1] - ylim[0] > 400.0:
            plot.ylim(0.0, None)

        if ref_plot_func is not None:
            ref_lines += ref_plot_func(color=profile_line[0].get_color())

        if len(ref_lines) > 0:
            plot.legend(ref_lines, [l.get_label() for l in ref_lines])

    plot.figlegend(
        lines,
        [l.get_label() for l in lines],
        loc="upper center",
        ncol=labels_ncol,
        numpoints=3,
        bbox_to_anchor=[0.5, 0.0],
        bbox_transform=plot.gcf().transFigure,
    )
    plot.grid(True)

    # remove height labels for plots not in first coloumn
    [fig.axes[n].set_ylabel("") for n in range(len(fig.axes)) if n % c != 0]
    [fig.axes[n].yaxis.set_ticklabels([]) for n in range(len(fig.axes)) if n % c != 0]

    title = "Vertical cloud profiles"
    if initial_condition is not None:
        title += "\n" + Var.repr(initial_condition)

    if all(
        [
            p.cloud_model.environment == profiles[0].cloud_model.environment
            for p in profiles
        ]
    ):
        title += "\nIn {}".format(str(profiles[0].cloud_model.environment))

    if all([p.cloud_model == profiles[0].cloud_model for p in profiles]):
        title += "\n{}".format(str(profiles[0].cloud_model))

    plot.suptitle(title, fontsize=14, y=1.1)

    return fig


def plot_hydrometeor_evolution(
    evolutions,
    variables=[
        "q_v",
    ],
    legend_loc="lower right",
    initial_condition=None,
    fig=None,
    ny=1,
    global_legend=False,
    labels=None,
    **kwargs
):
    n = len(variables)

    shape = int(math.ceil(float(n) / ny)), ny
    gs = gridspec.GridSpec(*shape)

    if fig is None:
        fig = plot.figure(figsize=(12, 6 * n))

    colors = {}
    for n, (v, s) in enumerate(zip(variables, list(gs))):
        lines = []

        t_max = 0.0
        data_max = 0.0
        data_min = 1.0e16
        plot.subplot(s)

        for n_evolution, evolution in enumerate(evolutions):
            try:
                i = Var.names.index(v)
                evolution_data = evolution.F[:, i]
                evolution_time = evolution.t
            except ValueError:
                if v in ["Nc", "r_c", "r_c__0"]:
                    try:
                        evolution_data = np.array(evolution.extra_vars[v])
                        evolution_time = evolution.extra_vars["t_substeps"]

                        if v in ["r_c", "r_c__0"]:
                            evolution_data *= 1.0e6
                    except KeyError:
                        continue
                else:
                    raise Exception("Variable %s not found" % v)

            if v in ["q_v", "q_l", "q_r"]:
                evolution_data = 1000.0 * np.array(evolution_data)

            if "markers" in kwargs:
                markers = kwargs["markers"]
                marker = markers[n_evolution]
            else:
                marker = "."

            if "colors" in kwargs:
                color = kwargs["colors"][n_evolution]
            else:
                color = colors.get(evolution, None)

            if not color is None:
                line = plot.plot(
                    evolution_time,
                    evolution_data,
                    label=str(evolution),
                    marker=marker,
                    linestyle="",
                    color=color,
                )
            else:
                line = plot.plot(
                    evolution_time,
                    evolution_data,
                    label=str(evolution),
                    marker=marker,
                    linestyle="",
                )
            if not v in colors:
                colors[evolution] = line[0].get_color()
            lines += line

            data_max = max(max(evolution_data), data_max)
            data_min = min(min(evolution_data), data_min)
            t_max = max(max(evolution.t), t_max)

            if v == "q_l":
                plot.ylabel("cloud water specific mass [g/kg]")
            if v == "q_v":
                plot.ylabel("water vapour specific mass [g/kg]")

                T = evolution.F[:, Var.T]
                p = evolution.F[:, Var.p]
                if hasattr(evolution.model, "qv_sat"):
                    q_v__sat = evolution.model.qv_sat(T=T, p=p)
                elif hasattr(evolution.model, "constants"):
                    model_parameterisations = (
                        parameterisations.ParametersationsWithSpecificConstants(
                            constants=evolution.model.constants
                        )
                    )
                    q_v__sat = model_parameterisations.pv_sat.qv_sat(T=T, p=p)
                else:
                    warnings.warn(
                        "Using default constant for calculating `qv_sat` for %s"
                        % str(evolution)
                    )
                    q_v__sat = parameterisations.pv_sat.qv_sat(T=T, p=p)
                q_v__sat = 1000.0 * np.array(q_v__sat)
                data_min = min(min(q_v__sat), data_min)
                data_max = max(max(q_v__sat), data_max)
                color = lines[n_evolution].get_color()
                plot.plot(
                    evolution.t,
                    q_v__sat,
                    marker="",
                    color=color,
                    label="",
                    linestyle=":",
                )
            elif v == "r_c" or v == "r_c__0":
                plot.ylabel("Cloud droplet radius [$\mu m$]")
            elif v == "Nc":
                plot.ylabel(r"Cloud droplet number concention [$m^{-3}$]")
            elif v == "T":
                plot.ylabel("Temperature [K]")
            elif v == "q_r":
                plot.ylabel("rain specific mass [g/kg]")
            elif v == "p":
                plot.ylabel("pressure [Pa]")

        d_data = 0.1 * (data_max - data_min)
        d_t = 0.1 * t_max

        plot.xlim(-d_t, t_max + d_t)
        plot.ylim(data_min - d_data, data_max + d_data)
        plot.xlabel("time [s]")
        plot.grid(True)

        if not global_legend:
            leg = plot.legend(loc=legend_loc)
            try:
                leg.get_frame().set_alpha(0.8)
            except AttributeError:
                pass

    title = "Hydrometeor evolution"
    if initial_condition is not None:
        title += "\n" + Var.repr(
            initial_condition,
            skip=[
                "r",
                "w",
                "z",
            ],
        )

    if global_legend:
        if labels is None:
            labels = [l.get_label() for l in lines]
        plot.figlegend(lines, labels, loc="lower center", ncol=ny, labelspacing=0.0)

    plot.suptitle(title)
    fig.subplots_adjust(top=0.95)

    return fig
