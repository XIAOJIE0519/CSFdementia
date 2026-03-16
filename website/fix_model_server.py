#!/usr/bin/env python3
import joblib, warnings
warnings.filterwarnings('ignore')

model = joblib.load('/www/wwwroot/csfdementia.top/data/machine/unified_LogisticRegression.pkl')
print('Model type:', type(model))

# Deep inspect all nested estimators
def fix_multi_class(obj, depth=0):
    fixed = 0
    if hasattr(obj, '__dict__'):
        if 'multi_class' in obj.__dict__:
            print(' ' * depth + f'Removing multi_class from {type(obj).__name__}')
            del obj.__dict__['multi_class']
            fixed += 1
    for attr in ['estimators_', 'estimator', 'base_estimator']:
        child = getattr(obj, attr, None)
        if child is None:
            continue
        if isinstance(child, list):
            for c in child:
                fixed += fix_multi_class(c, depth+2)
        else:
            fixed += fix_multi_class(child, depth+2)
    return fixed

n = fix_multi_class(model)
print(f'Fixed {n} multi_class attributes')

joblib.dump(model, '/www/wwwroot/csfdementia.top/data/machine/unified_LogisticRegression.pkl')
print('Saved.')

# Verify
import numpy as np
model2 = joblib.load('/www/wwwroot/csfdementia.top/data/machine/unified_LogisticRegression.pkl')
prob = model2.predict_proba(np.zeros((1,15)))[0][1]
print(f'Test OK: prob={prob:.4f}')
