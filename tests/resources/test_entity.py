import unittest

from mock import patch

from resources.entity import Entity


class TestEntity(unittest.TestCase):
    def setUp(self):
        return

    @patch('resources.entity.entities', [{'id': 1}])
    def test_get_200(self):
        # arrange
        test_input = 1
        expected = ({'entity': {'id': 1}}, 200)

        entity = Entity()

        # act
        actual = entity.get(test_input)

        # assert
        self.assertEqual(actual, expected)

    def test_get_404(self):
        # arrange
        test_input = 1
        expected = ({'entity': None}, 404)

        entity = Entity()

        # act
        actual = entity.get(test_input)

        # assert
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()