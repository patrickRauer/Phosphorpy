from datetime import datetime


class Log:

    _items = []

    def __init__(self):
        pass

    def add_entry(self, function_name, parameters):
        self._items.append(
            (
                datetime.now().isoformat(), function_name, parameters
            )
        )

    def write(self, path):
        with open(path, 'w') as f:
            for item in self._items:
                f.write('{}:\t {}\t{}\n'.format(item[0], item[1], item[2]))
