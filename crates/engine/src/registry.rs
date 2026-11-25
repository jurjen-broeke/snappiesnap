//! Command specifications and the registry that holds them.

use schemars::JsonSchema;
use serde::Serialize;
use serde::de::DeserializeOwned;
use serde_json::Value;

use crate::{CmdError, Session};

/// What a command may change. This decides how [`Session::execute`] protects and records it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Effect {
    /// Reads only. Gets `&Session`, so it can't change anything.
    Query,
    /// Changes the open molecule. Recorded for undo; restored on failure.
    Edit,
    /// Changes session state (open/close a document, undo/redo, selection). Not itself
    /// undoable; the document is restored on failure.
    Session,
}

type QueryFn = Box<dyn Fn(&Session, Value) -> Result<Value, CmdError> + Send + Sync>;
type MutateFn = Box<dyn Fn(&mut Session, Value) -> Result<Value, CmdError> + Send + Sync>;

pub(crate) enum Run {
    Query(QueryFn),
    Mutate(MutateFn),
}

/// Parameters for commands that take none. Accepts `{}` or `null`.
#[derive(Debug, Default, serde::Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NoParams {}

/// One command: its identity, where it appears, its parameter schema and what it does.
pub struct CommandSpec {
    /// Stable, `snake_case`; also the MCP tool name. Renaming one is a breaking change.
    pub id: &'static str,
    pub label: &'static str,
    pub description: &'static str,
    /// Menu path, for example `"Edit/Reverse Complement"`.
    pub menu: Option<&'static str>,
    /// Default shortcut, for example `"Cmd+Shift+R"`.
    pub shortcut: Option<&'static str>,
    pub effect: Effect,
    pub needs_document: bool,
    params_schema: Value,
    /// Valid example params. Every command needs at least one: they document the command for
    /// agents and seed the hostile-input suite, which mutates them field by field.
    pub examples: Vec<Value>,
    pub(crate) run: Run,
}

impl CommandSpec {
    /// A read-only command.
    pub fn query<P: DeserializeOwned + JsonSchema + 'static>(
        id: &'static str,
        label: &'static str,
        description: &'static str,
        run: fn(&Session, P) -> Result<Value, CmdError>,
    ) -> Self {
        Self::build::<P>(
            id,
            label,
            description,
            Effect::Query,
            true,
            Run::Query(Box::new(move |s, v| run(s, parse(id, v)?))),
        )
    }

    /// A command that edits the open molecule (undoable).
    pub fn edit<P: DeserializeOwned + JsonSchema + 'static>(
        id: &'static str,
        label: &'static str,
        description: &'static str,
        run: fn(&mut Session, P) -> Result<Value, CmdError>,
    ) -> Self {
        Self::build::<P>(
            id,
            label,
            description,
            Effect::Edit,
            true,
            Run::Mutate(Box::new(move |s, v| run(s, parse(id, v)?))),
        )
    }

    /// A command that changes session state.
    pub fn session<P: DeserializeOwned + JsonSchema + 'static>(
        id: &'static str,
        label: &'static str,
        description: &'static str,
        run: fn(&mut Session, P) -> Result<Value, CmdError>,
    ) -> Self {
        Self::build::<P>(
            id,
            label,
            description,
            Effect::Session,
            false,
            Run::Mutate(Box::new(move |s, v| run(s, parse(id, v)?))),
        )
    }

    fn build<P: JsonSchema>(
        id: &'static str,
        label: &'static str,
        description: &'static str,
        effect: Effect,
        needs_document: bool,
        run: Run,
    ) -> Self {
        CommandSpec {
            id,
            label,
            description,
            menu: None,
            shortcut: None,
            effect,
            needs_document,
            params_schema: serde_json::to_value(schemars::schema_for!(P)).unwrap_or(Value::Null),
            examples: Vec::new(),
            run,
        }
    }

    /// Add a valid example of the params.
    pub fn example(mut self, params: Value) -> Self {
        self.examples.push(params);
        self
    }

    pub fn menu(mut self, path: &'static str) -> Self {
        self.menu = Some(path);
        self
    }

    pub fn shortcut(mut self, keys: &'static str) -> Self {
        self.shortcut = Some(keys);
        self
    }

    pub fn needs_document(mut self, needs: bool) -> Self {
        self.needs_document = needs;
        self
    }

    /// JSON Schema of the params object.
    pub fn params_schema(&self) -> &Value {
        &self.params_schema
    }

    /// Everything about the command except how it runs.
    pub fn info(&self, with_schema: bool) -> CommandInfo {
        CommandInfo {
            id: self.id,
            label: self.label,
            description: self.description,
            menu: self.menu,
            shortcut: self.shortcut,
            effect: self.effect,
            needs_document: self.needs_document,
            params: with_schema.then(|| self.params_schema.clone()),
            examples: if with_schema {
                self.examples.clone()
            } else {
                Vec::new()
            },
        }
    }
}

/// The serialisable description of a command (`commands.list`, `commands.describe`).
#[derive(Debug, Clone, Serialize)]
pub struct CommandInfo {
    pub id: &'static str,
    pub label: &'static str,
    pub description: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub menu: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shortcut: Option<&'static str>,
    pub effect: Effect,
    pub needs_document: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub params: Option<Value>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub examples: Vec<Value>,
}

fn parse<P: DeserializeOwned>(id: &str, params: Value) -> Result<P, CmdError> {
    let params = if params.is_null() {
        Value::Object(serde_json::Map::new())
    } else {
        params
    };
    serde_json::from_value(params).map_err(|e| {
        CmdError::invalid_params(format!("{id}: {e}")).with_hint(format!(
            "`commands.describe` with id \"{id}\" shows the parameters"
        ))
    })
}

/// All commands, by id.
pub struct Registry {
    commands: Vec<CommandSpec>,
}

impl Registry {
    /// Build a registry. Commands with a malformed or duplicate id are left out and reported.
    pub fn new(specs: Vec<CommandSpec>) -> (Registry, Vec<String>) {
        let mut problems = Vec::new();
        let mut commands: Vec<CommandSpec> = Vec::with_capacity(specs.len());
        for spec in specs {
            if !valid_id(spec.id) {
                problems.push(format!(
                    "command id {:?} must be 1-64 chars of a-z, 0-9 and _",
                    spec.id
                ));
            } else if commands.iter().any(|c| c.id == spec.id) {
                problems.push(format!("duplicate command id {:?}", spec.id));
            } else {
                commands.push(spec);
            }
        }
        (Registry { commands }, problems)
    }

    pub fn get(&self, id: &str) -> Option<&CommandSpec> {
        self.commands.iter().find(|c| c.id == id)
    }

    pub fn iter(&self) -> impl Iterator<Item = &CommandSpec> {
        self.commands.iter()
    }

    pub fn len(&self) -> usize {
        self.commands.len()
    }

    pub fn is_empty(&self) -> bool {
        self.commands.is_empty()
    }
}

fn valid_id(id: &str) -> bool {
    (1..=64).contains(&id.len())
        && id
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'_')
        && id.bytes().next().is_some_and(|b| b.is_ascii_lowercase())
}
