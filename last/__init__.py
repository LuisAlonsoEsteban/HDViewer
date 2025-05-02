
PARAMS = {
    'color': {
        'nan': [0.5, 0.5, 0.5, 0.4]
    },
    'edgecolor': {
        'nan': "grey"
    },
    'linewidth': {
        'min': 0.5,
        'max': 5.5,
        'nan': 0.2
    },
    'alpha': {
        'min': 0.5,
        'max': 1.0,
        'nan': 0.1
    }
}

RAM_LIMIT = 0.1  # GB

class Config:
    def __init__(self, defaults):
        self.settings = defaults.copy()

    def set(self, key, value):
        self.settings[key] = value

    def get(self, key):
        return self.settings.get(key)

    def update(self, new_configs):
        self.settings.update(new_configs)

    def reset(self):
        self.settings = PARAMS.copy()

    def __getitem__(self, key):
        return self.settings[key]
    
    def __setitem__(self, key, value):
        self.settings[key] = value

config = Config(PARAMS)
