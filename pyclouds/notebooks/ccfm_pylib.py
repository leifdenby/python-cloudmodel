def moist_adjust(tem, q_v, prs):
    q_sat = lua( max( 180._dp, min( tem, 370._dp))) / prs
    q_sat = min( 0.5_dp, q_sat )
    cor   = 1._dp / ( 1._dp - ( Rv/Rd - 1._dp ) * q_sat )
    q_sat = q_sat * cor

    cond  = (q_v - q_sat) (1._dp + q_sat * cor * lub( max(180._dp,min( tem ,370._dp))))
    cond  =  max( cond, 0._dp )
    tem = tem + luc(tem) * cond
    q_v  = q_v - cond

    if ( cond > 0._dp ):
        q_sat = lua( max( 180._dp, min( tem, 370._dp))) / prs
        q_sat = min( 0.5_dp, q_sat )
        cor   = 1._dp / ( 1._dp - ( Rv/Rd - 1._dp ) * q_sat )
        q_sat = q_sat * cor
        cond  = (q_v - q_sat) /                                             &
        (1._dp + q_sat * cor * lub( max(180._dp,min(tem,370._dp))) )
        cond  =  max( cond, 0._dp )
        tem = tem + luc(tem) * cond
        q_v  = q_v - cond

    return q_v


def lua(:q
        )
