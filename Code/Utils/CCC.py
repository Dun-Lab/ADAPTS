def filter_func_cellphone_pvals(df, key='cellphone_pvals', thres=0.05):
    return df[key] <= thres

def filter_func_specificity_rank(df, key='specificity_rank', thres=0.05):
    return df[key] <= thres
