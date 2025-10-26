# Developer Best Practices

- [ ] I remember that code quality in ADA *IS GRADED!*
- [ ] I don't use `logger.info()`; instead I always use the logger from `src/utils/logging_utils.py`, for example:

    ```python
    from src.utils.logging_utils import logger
    logger.info("Hello, world!")
    ```

- [ ] I always install (`pre-commit install`) pre-commit hooks and never skip them when committing.
- [ ] I always leave plenty of comments in the code.
- [ ] Every public class, method, file, and function has a docstring.
- [ ] I always use type annotations in function signatures, for example:

    ```python
    def add(a: int, b: int) -> int:
        return a + b
    ```

- [ ] Whenever I commit code, I also add unit and integration tests for any new functionality.
- [ ] I never run `git add .`; instead, I always run `git status` to see what files have changed and then add them one by one.
- [ ] I never commit secrets, API keys, or irrelevant files such as duplicates, caches, boilerplate, etc.
